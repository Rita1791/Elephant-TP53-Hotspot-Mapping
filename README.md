from pathlib import Path
import zipfile, textwrap

root = Path("/mnt/data/Elephant-TP53-GitHub-Research-Update")
(root / "docs").mkdir(parents=True, exist_ok=True)

readme = r'''<div align="center">

# 🧬 Comparative In-Silico Analysis of TP53 Mutation Hotspots Between Humans and Elephants

### Exploring evolutionary conservation, retrogene variation, and TP53 mutation hotspots in the context of Peto's paradox

[![Research](https://img.shields.io/badge/Research-Comparative%20Genomics-6A1B9A?style=for-the-badge)](https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping)
[![Cancer Biology](https://img.shields.io/badge/Cancer%20Biology-TP53-B71C1C?style=for-the-badge)](https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Computational%20Biology-00599C?style=for-the-badge)](https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping)
[![Status](https://img.shields.io/badge/Status-Completed-2E7D32?style=for-the-badge)](https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping)
[![Publication](https://img.shields.io/badge/Publication-March%202026-0A7EA4?style=for-the-badge)](https://doi.org/10.25215/9141002199)

</div>

---

## 🔎 Research at a Glance

| | |
|---|---|
| **Research area** | Comparative genomics • Cancer biology • Evolutionary bioinformatics |
| **Primary gene** | **TP53 / p53** |
| **Comparative system** | Human • African elephant (*Loxodonta africana*) • Asian elephant (*Elephas maximus*) |
| **Focal sites** | R175 • G245 • R248 • R249 • R273 • R282 |
| **Approach** | Protein-level multiple sequence alignment • hotspot mapping • conservation analysis |
| **Key tools** | NCBI • UniProt • MAFFT • Jalview • MEGA • Python • Google Colab |
| **Research status** | ✅ Completed |
| **Publication** | RED'SHINE Publication, March 2026 |

---

## 📄 Abstract

This research presents a comparative in-silico analysis of major human **TP53 mutation hotspots** in elephant canonical and retrogene sequences. The study was motivated by **Peto's paradox**: despite their large body size and long lifespan, elephants do not experience cancer incidence proportional to the number of cells or years at risk.

Six clinically significant human TP53 mutation hotspots — **R175, G245, R248, R249, R273, and R282** — were mapped onto African and Asian elephant TP53 protein sequences using protein-level multiple sequence alignment. Publicly available sequence resources from **NCBI** and **UniProt** were analysed using **MAFFT, Jalview, MEGA, and Python-based workflows implemented in Google Colab**.

The analysis indicates strong conservation of key DNA-binding residues in canonical elephant TP53 sequences, whereas retrogene copies exhibit position-specific variability. A descriptive composite scoring framework was additionally used to prioritise hotspots using mutation frequency, cross-species conservation, and retrogene variability.

The project demonstrates how transparent comparative genomics workflows can be used to investigate evolutionary constraints surrounding tumour-suppressor genes and generate hypotheses for future experimental investigation.

---

## 🧠 Biological Rationale

**TP53** encodes the tumour-suppressor protein p53, a central regulator of genomic stability, DNA-damage response, cell-cycle arrest, senescence, and apoptosis.

Human cancers contain recurrent TP53 mutations, particularly within functionally important regions of the protein. These recurrent hotspots provide a useful reference for asking whether residues frequently altered in human cancer are conserved across species.

Elephants provide an important comparative system because their large body size and long lifespan create a biological context in which the expected relationship between cellular burden and cancer risk is not straightforward. This phenomenon is widely discussed as **Peto's paradox**.

Rather than assuming that TP53 sequence conservation directly explains cancer resistance, this study uses comparative sequence analysis to identify patterns of conservation and divergence that can inform future functional research.

---

## 🎯 Aim & Objectives

### Aim

To perform a structured in-silico comparative mapping of major human TP53 mutation hotspots across canonical and retrogene TP53 protein sequences in elephants.

### Objectives

1. Retrieve curated TP53 protein sequences from humans, African elephants, and Asian elephants using publicly accessible biological databases.
2. Perform protein-level multiple sequence alignment to establish positional correspondence between sequences.
3. Identify and map major human TP53 mutation hotspots: **R175, G245, R248, R249, R273, and R282**.
4. Evaluate conservation and amino-acid variation at corresponding positions.
5. Compare canonical elephant TP53 sequences with TP53 retrogene-derived sequences.
6. Apply a transparent descriptive scoring framework integrating mutation frequency, conservation, and retrogene variability.
7. Establish a documented computational workflow suitable for inspection and future extension.

---

## 🗂️ Data Sources

| Resource | Role in the study |
|---|---|
| 🧬 **NCBI** | Sequence retrieval and biological database resources |
| 🧬 **UniProt** | Curated protein sequence and annotation resources |
| 🧬 **Human TP53 mutation data** | Identification of recurrent cancer-associated hotspots |
| 🐘 **Elephant TP53 sequences** | Cross-species sequence comparison |
| 🧬 **Elephant TP53 retrogenes** | Assessment of sequence-level variability |

No patient-identifiable genomic data are required for this comparative study.

---

## 🔬 Computational Methodology

```text
Human TP53 Mutation Hotspots
              │
              ▼
     Sequence Retrieval
              │
              ▼
       Data Preprocessing
              │
              ▼
 Human–Elephant TP53 Comparison
              │
              ▼
 Multiple Sequence Alignment (MAFFT)
              │
              ▼
      Hotspot Residue Mapping
              │
              ▼
 Conservation & Variation Analysis
              │
              ▼
 Canonical vs Retrogene Comparison
              │
              ▼
     Composite Descriptive Score
              │
              ▼
     Biological Interpretation
