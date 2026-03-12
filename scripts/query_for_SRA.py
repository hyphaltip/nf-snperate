#!/usr/bin/env python3
"""
query_for_SRA.py

Query the NCBI BioSample (or SRA) database for sequencing runs, apply
metadata filters, and write results to an incrementally-updated CSV so
that periodic re-runs only append new accessions.

Requirements:
    pip install biopython

Basic usage:
    python query_for_SRA.py \\
        --query "Aspergillus fumigatus" \\
        --sra_filter "DNA" --sra_filter "paired" \\
        --email you@example.com

With date filtering and API key:
    python query_for_SRA.py \\
        --query "Aspergillus fumigatus" \\
        --sra_filter "WGS" --sra_filter "ILLUMINA" \\
        --since 2023-01-01 \\
        --output sra_results.csv \\
        --email you@example.com \\
        --api_key YOUR_NCBI_API_KEY

Filter terms are matched case-insensitively against these run metadata
fields: LibrarySource, LibraryStrategy, LibraryLayout, LibrarySelection,
Platform, Model.

Common shortcuts (expanded automatically):
    DNA       → matches "GENOMIC" in LibrarySource
    RNA       → matches "TRANSCRIPTOMIC" in LibrarySource
    genome    → matches "WGS" in LibraryStrategy
    paired    → matches "PAIRED" in LibraryLayout
    single    → matches "SINGLE" in LibraryLayout
    Illumina  → matches "ILLUMINA" in Platform
    ONT       → matches "OXFORD_NANOPORE" in Platform
    PacBio    → matches "PACBIO_SMRT" in Platform
"""

import argparse
import csv
import io
import logging
import os
import sys
import time
from datetime import datetime
from typing import Dict, List, Optional, Set, Tuple

from Bio import Entrez

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Batch size for efetch calls (NCBI recommends ≤ 500)
BATCH_SIZE = 200

# Minimum sleep between requests (NCBI rate limits)
#   Without API key: 3 req/s → 0.34 s
#   With API key:   10 req/s → 0.11 s
SLEEP_NO_KEY  = 0.34
SLEEP_API_KEY = 0.11

# Preferred column ordering in the output CSV
OUTPUT_COLUMNS = [
    "Run",
    "BioSample",
    "BioProject",
    "SRAStudy",
    "Experiment",
    "ScientificName",
    "TaxID",
    "SampleName",
    "LibraryName",
    "LibraryStrategy",
    "LibrarySelection",
    "LibrarySource",
    "LibraryLayout",
    "Platform",
    "Model",
    "spots",
    "bases",
    "size_MB",
    "avgLength",
    "spots_with_mates",
    "CenterName",
    "Submission",
    "ReleaseDate",
    "LoadDate",
]

# Synonym expansion: user-supplied term → list of substrings to test against
# the combined metadata string (case-insensitive).
_FILTER_SYNONYMS: Dict[str, List[str]] = {
    "dna":      ["genomic"],
    "rna":      ["transcriptomic"],
    "genome":   ["wgs"],
    "wgs":      ["wgs"],
    "paired":   ["paired"],
    "single":   ["single"],
    "illumina": ["illumina"],
    "ont":      ["oxford_nanopore"],
    "nanopore": ["oxford_nanopore"],
    "pacbio":   ["pacbio_smrt"],
    "amplicon": ["amplicon"],
    "rnaseq":   ["rna-seq"],
    "chip":     ["chip-seq"],
    "atac":     ["atac-seq"],
}

# Mapping from user filter term to NCBI esearch field tag, used to narrow
# the server-side query before client-side filtering.
_NCBI_FILTER_TAGS: Dict[str, str] = {
    "dna":      '"GENOMIC"[Source]',
    "rna":      '"TRANSCRIPTOMIC"[Source]',
    "genome":   '"wgs"[Strategy]',
    "wgs":      '"wgs"[Strategy]',
    "paired":   '"paired"[Properties]',
    "single":   '"single"[Properties]',
    "illumina": '"illumina"[Platform]',
    "ont":      '"OXFORD_NANOPORE"[Platform]',
    "nanopore": '"OXFORD_NANOPORE"[Platform]',
    "pacbio":   '"PACBIO_SMRT"[Platform]',
    "amplicon": '"amplicon"[Strategy]',
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def expand_filter_term(term: str) -> List[str]:
    """Return the list of lowercase substrings that represent *term*."""
    return _FILTER_SYNONYMS.get(term.lower(), [term.lower()])


def run_matches_filters(run: Dict, filters: List[str]) -> bool:
    """
    Return True when *run* satisfies **all** filter terms (AND logic).

    Each term is tested as a substring against the combined lowercase text of
    LibrarySource, LibraryStrategy, LibraryLayout, LibrarySelection, Platform,
    and Model.
    """
    if not filters:
        return True
    haystack = " ".join(
        run.get(col, "") for col in
        ("LibrarySource", "LibraryStrategy", "LibraryLayout",
         "LibrarySelection", "Platform", "Model")
    ).lower()
    for term in filters:
        alternatives = expand_filter_term(term)
        if not any(alt in haystack for alt in alternatives):
            return False
    return True


def build_ncbi_query(dbname: str, query: str, sra_filters: List[str], since: Optional[str]) -> str:
    """
    Compose a full NCBI esearch query string by combining the user query,
    any filter terms that have known NCBI field tags, and an optional date
    range.  Unknown filter terms are left for client-side filtering only.
    """
    parts = [query]
    seen_tags: Set[str] = set()
    if dbname == "sra":
        for f in sra_filters:        
            tag = _NCBI_FILTER_TAGS.get(f.lower())
            print(f"Filter term: {f}  →  NCBI tag: {tag}")
            if tag and tag not in seen_tags:
                parts.append(tag)
                seen_tags.add(tag)
    if since:
        since_fmt = since.replace("-", "/")
        parts.append(f'("{since_fmt}"[PDAT] : "3000/01/01"[PDAT])')
    return " AND ".join(parts)


def load_existing(output_path: str) -> Tuple[Set[str], List[str]]:
    """
    Read the output CSV (if it exists) and return:
        (set of run accessions already present, list of column headers)
    """
    if not os.path.exists(output_path):
        return set(), []
    accessions: Set[str] = set()
    fieldnames: List[str] = []
    with open(output_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            acc = (row.get("Run") or "").strip()
            if acc:
                accessions.add(acc)
    return accessions, fieldnames


def parse_runinfo(raw) -> List[Dict]:
    """Parse the SRA runinfo CSV text returned by efetch into a list of dicts.
    Accepts both str and bytes (some Biopython versions return bytes)."""
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8", errors="replace")
    rows = []
    try:
        reader = csv.DictReader(io.StringIO(raw))
        for row in reader:
            run = row.get("Run", "").strip()
            if run and not run.startswith("<"):  # guard against XML error blobs
                rows.append(dict(row))
    except Exception as exc:
        logging.warning("Could not parse runinfo: %s", exc)
        logging.debug("Response (first 500 chars): %s", raw[:500])
    return rows


def safe_efetch(sleep: float, **kwargs) -> str:
    """Call Entrez.efetch with up to 3 retries; always returns str."""
    for attempt in range(3):
        try:
            handle  = Entrez.efetch(**kwargs)
            content = handle.read()
            handle.close()
            if isinstance(content, bytes):
                content = content.decode("utf-8", errors="replace")
            time.sleep(sleep)
            return content
        except Exception as exc:
            wait = sleep * (4 ** attempt)
            logging.warning("efetch attempt %d failed: %s — retrying in %.1fs",
                            attempt + 1, exc, wait)
            time.sleep(wait)
    raise RuntimeError("efetch failed after 3 attempts")


def safe_elink_batch(ids: List[str], dbfrom: str, db: str,
                     sleep: float) -> List[str]:
    """
    Run elink for a batch of BioSample IDs and return the linked SRA IDs.
    Returns an empty list on failure (logged as a warning).
    """
    id_str = ",".join(ids)
    for attempt in range(3):
        try:
            print(f"Linking BioSample IDs: {id_str}")
            handle = Entrez.elink(dbfrom=dbfrom, db=db, id=id_str)
            records = Entrez.read(handle)
            handle.close()
            time.sleep(sleep)
            linked: List[str] = []
            for record in records:
                for link_set in record.get("LinkSetDb", []):
                    linked.extend(str(l["Id"]) for l in link_set.get("Link", []))
            return linked
        except Exception as exc:
            wait = sleep * (10 ** attempt)
            logging.warning("elink attempt %d failed: %s — retrying in %.1fs",
                            attempt + 1, exc, wait)
            time.sleep(wait)
    logging.error("elink failed for batch of %d IDs after 3 attempts", len(ids))
    return []


# ---------------------------------------------------------------------------
# Search strategies
# ---------------------------------------------------------------------------

def search_via_sra(query_str: str, filters: List[str],
                   max_results: int, sleep: float) -> List[Dict]:
    """
    Search SRA directly, fetch runinfo, apply client-side filters.
    Returns a list of run-info dicts.
    """
    logging.info("Searching SRA with query: %s", query_str)
    handle = Entrez.esearch(db="sra", term=query_str,
                            usehistory="y", retmax=0)
    rec = Entrez.read(handle)
    handle.close()
    time.sleep(sleep)

    total = int(rec["Count"])
    webenv    = rec["WebEnv"]
    query_key = rec["QueryKey"]
    logging.info("SRA hits: %d", total)
    if total == 0:
        return []
    if max_results:
        total = min(total, max_results)

    runs: List[Dict] = []
    for start in range(0, total, BATCH_SIZE):
        batch = min(BATCH_SIZE, total - start)
        logging.info("  fetching %d-%d / %d ...", start + 1, start + batch, total)
        content = safe_efetch(
            sleep,
            db="sra", rettype="runinfo", retmode="text",
            retstart=start, retmax=batch,
            webenv=webenv, query_key=query_key,
        )
        runs.extend(parse_runinfo(content))

    logging.info("Fetched %d run records from SRA", len(runs))
    return runs


def search_via_biosample(query_str: str, filters: List[str],
                        max_results: int, sleep: float) -> List[Dict]:
    """
    Search BioSample, link to SRA, fetch runinfo, apply client-side filters.
    Returns a list of run-info dicts.
    """
    logging.info("Searching BioSample with query: %s", query_str)
    handle = Entrez.esearch(db="biosample", term=query_str,
                            usehistory="y", retmax=0)
    rec = Entrez.read(handle)
    handle.close()
    time.sleep(sleep)

    total_bs = int(rec["Count"])
    webenv    = rec["WebEnv"]
    query_key = rec["QueryKey"]
    logging.info("BioSample hits: %d", total_bs)
    if total_bs == 0:
        return []
    if max_results:
        total_bs = min(total_bs, max_results)

    # Link BioSample search history directly to SRA using webenv/query_key.
    # This avoids fetching BioSample UIDs as an intermediate step (which was
    # the source of the "#exLinkSrv2 address table is empty" elink error).
    # linkname="biosample_sra" pins the exact NCBI link type and prevents the
    # link-server resolution failure that occurs when linkname is left unset.
    logging.info("Linking BioSample results to SRA (linkname=biosample_sra) ...")
    for attempt in range(3):
        try:
            link_handle = Entrez.elink(
                dbfrom="biosample", db="sra",
                webenv=webenv, query_key=query_key,
                linkname="biosample_sra",
                cmd="neighbor",
            )
            link_recs = Entrez.read(link_handle)
            link_handle.close()
            time.sleep(sleep)
            break
        except Exception as exc:
            wait = sleep * (4 ** attempt)
            logging.warning("elink attempt %d failed: %s — retrying in %.1fs",
                            attempt + 1, exc, wait)
            time.sleep(wait)
    else:
        logging.error("elink BioSample->SRA failed after 3 attempts")
        return []

    sra_ids: List[str] = []
    for record in link_recs:
        for link_set in record.get("LinkSetDb", []):
            for link in link_set.get("Link", []):
                sra_ids.append(str(link["Id"]))

    sra_ids = list(dict.fromkeys(sra_ids))  # deduplicate, preserve order
    logging.info("Linked to %d SRA experiment IDs", len(sra_ids))
    if not sra_ids:
        return []

    # Fetch runinfo for each SRA ID batch
    runs: List[Dict] = []
    for i in range(0, len(sra_ids), BATCH_SIZE):
        batch = sra_ids[i : i + BATCH_SIZE]
        logging.info("  fetching runinfo for SRA IDs %d-%d / %d ...",
                     i + 1, i + len(batch), len(sra_ids))
        content = safe_efetch(
            sleep,
            db="sra", rettype="runinfo", retmode="text",
            id=",".join(batch),
        )
        runs.extend(parse_runinfo(content))

    logging.info("Fetched %d run records via BioSample link", len(runs))
    return runs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Query NCBI BioSample/SRA and save run accessions to CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Requirements:")[0].strip(),
    )
    parser.add_argument(
        "--query", required=True,
        help='Search term, e.g. "Aspergillus fumigatus[Organism]"',
    )
    parser.add_argument(
        "--sra_filter", action="append", default=[], dest="sra_filters",
        metavar="TERM",
        help="Metadata filter term (repeatable); all terms must match (AND). "
             "See module docstring for shortcuts.",
    )
    parser.add_argument(
        "--since", default=None, metavar="YYYY-MM-DD",
        help="Only include runs with a publication date on or after this date.",
    )
    parser.add_argument(
        "--output", default="sra_results.csv",
        help="Output CSV file (default: sra_results.csv). "
             "Existing rows are preserved; only new accessions are appended.",
    )
    parser.add_argument(
        "--email", required=True,
        help="E-mail address required by NCBI Entrez policy.",
    )
    parser.add_argument(
        "--api_key", default=None,
        help="NCBI API key for higher rate limits (10 req/s vs 3 req/s).",
    )
    parser.add_argument(
        "--db", choices=["biosample", "sra"], default="biosample",
        help="Primary database to search. "
             "'biosample' (default) searches BioSample then links to SRA. "
             "'sra' searches SRA directly (faster for very large queries).",
    )
    parser.add_argument(
        "--max_results", type=int, default=0, metavar="N",
        help="Stop after N results from NCBI (0 = no limit).",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Enable DEBUG-level logging.",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )

    # Validate date
    if args.since:
        try:
            datetime.strptime(args.since, "%Y-%m-%d")
        except ValueError:
            parser.error(f"--since must be YYYY-MM-DD, got: {args.since!r}")

    # Configure Entrez
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    sleep = SLEEP_API_KEY if args.api_key else SLEEP_NO_KEY

    # Load existing results for incremental update
    existing_accs, existing_cols = load_existing(args.output)
    logging.info("Accessions already in %s: %d", args.output, len(existing_accs))

    # Build NCBI query
    query_str = build_ncbi_query(args.db, args.query, args.sra_filters, args.since)
    logging.info("NCBI query: %s", query_str)

    # Fetch run records
    if args.db == "biosample":
        all_runs = search_via_biosample(
            query_str, args.sra_filters, args.max_results, sleep)
    else:
        all_runs = search_via_sra(
            query_str, args.sra_filters, args.max_results, sleep)

    if not all_runs:
        logging.info("No runs returned from NCBI.")
        sys.exit(0)

    # Apply client-side filters and skip already-seen accessions
    new_runs: List[Dict] = []
    skipped_existing  = 0
    skipped_filter    = 0
    for run in all_runs:
        acc = run.get("Run", "").strip()
        if not acc:
            continue
        if acc in existing_accs:
            skipped_existing += 1
            continue
        if not run_matches_filters(run, args.sra_filters):
            skipped_filter += 1
            logging.debug("Filtered out %s (%s / %s / %s)",
                          acc,
                          run.get("LibrarySource", "?"),
                          run.get("LibraryStrategy", "?"),
                          run.get("LibraryLayout", "?"))
            continue
        new_runs.append(run)
        existing_accs.add(acc)  # prevent duplicates within this run

    logging.info(
        "Result summary - total fetched: %d  |  already in file: %d  "
        "|  filtered out: %d  |  new: %d",
        len(all_runs), skipped_existing, skipped_filter, len(new_runs),
    )

    if not new_runs:
        logging.info("Nothing new to write.")
        sys.exit(0)

    # Determine column set
    # If appending, reuse the existing header; otherwise build from data + preferred order
    if existing_cols:
        fieldnames = existing_cols
        # Append any new columns the current fetch introduced
        for col in new_runs[0].keys():
            if col not in fieldnames:
                fieldnames.append(col)
    else:
        all_keys = list(new_runs[0].keys())
        preferred = [c for c in OUTPUT_COLUMNS if c in all_keys]
        extra     = [c for c in all_keys        if c not in preferred]
        fieldnames = preferred + extra

    # Append to output CSV (write header only when creating a new file)
    write_header = not os.path.exists(args.output)
    with open(args.output, "a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames,
                                extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerows(new_runs)

    logging.info("Appended %d new rows to %s", len(new_runs), args.output)


if __name__ == "__main__":
    main()
