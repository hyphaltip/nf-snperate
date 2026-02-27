# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**snperate** is a Nextflow-based SNP calling and population genomics pipeline for fungal pathogens (Aspergillus fumigatus). It performs alignment, variant calling, CNV detection, and phylogenetic analysis.

## Common Commands

```bash
# Run the Nextflow pipeline
nextflow run nextflow.nf -c config.txt

# Run individual pipeline steps (when using shell-based workflow)
# Index reference genome
./pipeline/00_index.sh
# Align reads
./pipeline/01_align.sh
# Call GVCF
./pipeline/02_call_gvcf.sh
# Joint genotyping + VCF filtering
./pipeline/03_jointGVCF_call_slice.sh
# or with bcftools:
./pipeline/03_jointGVCF_call_slice_bcftools.sh
# Combine VCFs
./pipeline/04_combine_vcf.sh
# CNV detection with mosdepth
./pipeline/11_mosdepth.sh
```

## Configuration

- `config.txt` - Main pipeline configuration (paths, tools, parameters)
- `samples.csv` - Sample manifest: `sample_id,fastq_paths(; delimited),platform`
- `population_sets.yaml` - Population definitions for subsetting VCFs by strain groups

## Architecture

### Workflow (Nextflow)
The main entry point is `nextflow.nf` which imports `SNPERATE_PIPELINE` from `./workflows/snperate`. The workflows directory needs to be created.

### Pipeline Steps
1. **Indexing** (00) - Build bwa index
2. **Alignment** (01) - bwa-mem2 alignment → CRAM files
3. **GVCF Calling** (02) - GATK HaplotypeCaller per-sample
4. **Joint Genotyping** (03) - Combine GVCFs and genotype:
   - GATK-based: `03_jointGVCF_call_slice.sh`
   - bcftools-based: `03_jointGVCF_call_slice_bcftools.sh`
5. **VCF Combination** (04) - Merge multiple VCFs
6. **Unmapped Assembly** (05) - Assemble unmapped reads
7. **Phylogenetics** (06-07) - SNP matrix → IQ-Tree
8. **Annotation** (08) - snpEff variant annotation
9. **CNV Detection** (11) - mosdepth coverage windows + R plotting

### Data Flow
```
samples.csv → FASTQ → bwa → CRAM → GATK → GVCF → JointGT → VCF → Filtered subsets by population
                                        ↓
                                  Phylogeny (IQ-Tree)
                                  CNV (mosdepth + R)
```

## Population Subsetting

VCF filtering by population uses `population_sets.yaml`. Each entry lists strain IDs belonging to a population. The pipeline can generate population-specific VCF subsets for downstream analysis.

## Development

- Use ruff and precommit for linting (per AGENTS.md)
- Python scripts in `scripts/` handle auxiliary tasks
- R scripts for CNV plotting and phylogenetics
