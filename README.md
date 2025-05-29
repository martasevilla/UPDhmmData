# UPDhmm Simulation and Evaluation Framework

This repository contains the complete pipeline for simulating and detecting Uniparental Disomy (UPD) events in trio-based next-generation sequencing (NGS) data. It supports both exome and genome simulations and includes comparisons with state-of-the-art methods (UPDio and AltAFplotter).

> **UPDhmm: Detecting Uniparental Disomy through NGS Trio Data**  
> Marta Sevilla Porras, Carlos Ruiz Arenas*, Luis Pérez Jurado*  
> *These authors equally contributed

**Bioconductor package**: [UPDhmm](https://www.bioconductor.org/packages/release/bioc/html/UPDhmm.html)  
**Source code**: [https://github.com/martasevilla/UPDhmm](https://github.com/martasevilla/UPDhmm)  

---

## Overview

Uniparental disomy (UPD) is a chromosomal anomaly where both homologs are inherited from the same parent. UPDhmm is an R/Bioconductor package that uses a Hidden Markov Model (HMM) to detect UPD events and classify them into isodisomy or heterodisomy. This repository provides:

- Scripts for generating simulated UPD events in real data (exome and genome).
- Scripts to run UPDhmm, UPDio, and AltAFplotter.
- Scripts to compute false positive rates using clean trio data.
- SSC cohort analysis described in the manuscript.

---

## Folder Structure

| Script | Description |
|--------|-------------|
| `01a-createDataExome.sh` | Prepares GIAB exome trio VCFs (region and SV filtering). |
| `01b-createDataGenome.sh` | Prepares high-coverage genome VCFs from 1000 Genomes Project. |
| `02-createRegions.sh` | Selects random genomic regions of specific sizes (1–20 Mb). |
| `03-createSimulations.sh` | Simulates 100 UPD events (5 per type/size) using real VCFs. |
| `04-runUPDio.sh` | Runs UPDio on simulated multi-sample VCFs. |
| `05-runUPDhmm.sh` | Runs UPDhmm and outputs predicted UPD regions. |
| `06-runAltAFplotter.sh` | Runs bcftools + Python summarization for AltAFplotter strategy. |
| `07-calculateFalsePositives.sh` | Applies all three methods on non-simulated data to compute false positive rates. |
| `08-processSFARI.sh` | Preprocesses WGS VCFs from the Simons Simplex Collection (SFARI) and calcualte UPD events using UPDhmm. |


---

## Requirements

- **R** ≥ 4.2 with Bioconductor packages:
  - `UPDhmm`, `VariantAnnotation`
- **bcftools** ≥ 1.16
- **Python** ≥ 3.9 with:
  - `pandas`, `cyvcf2`, `altair`, `openpyxl`, `multiprocess`
- **Perl** and the UPDio script

---

## How to Run

### Step 1: Data Preparation

```bash
# For exome data (GIAB)
01a-createDataExome.sh

# For genome data (1000 Genomes Project)
01ab-createDataGenome.sh
```

### Step 2: Region Selection

```bash
02-createRegions.sh
```

### Step 3: UPD Simulation

```bash
03-createSimulations.sh (sbatch --array=1-100)
```

### Step 4: Detection Methods

```bash
04-runUPDio.sh
05-runUPDhmm.sh
06-runAltAFplotter.sh
```

### Step 5: False Positive Estimation

```bash
07-calculateFalsePositives.sh <data/01_centro_seg_dup_SV/>
```
```bash
08-processSFARI.sh
```


---

## Reference Dataset

- **Exome**: GIAB trios (NA19675, NA19685)
- **Genome**: 10 high-coverage trios from 1000 Genomes
- **Real cohort**: 4307 trios from Simons Simplex Collection (SSC)

---

## Citation

If you use this codebase or UPDhmm in your research, please cite:

> Marta Sevilla Porras, Carlos Ruiz Arenas, Luis Pérez Jurado.  
> **UPDhmm: detecting Uniparental Disomy through NGS trio data**

---

## License

This project is distributed under the MIT License.
