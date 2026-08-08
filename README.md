<div align="center">

# 🧬 Comparative In-Silico Analysis of TP53 Mutation Hotspots Between Humans and Elephants

### Exploring evolutionary conservation, retrogene variation, and TP53 mutation hotspots in the context of Peto's paradox

<p align="center">

<img src="https://img.shields.io/badge/Research-Comparative%20Genomics-6A1B9A?style=for-the-badge">

<img src="https://img.shields.io/badge/Cancer%20Biology-TP53-B71C1C?style=for-the-badge">

<img src="https://img.shields.io/badge/Bioinformatics-Computational%20Biology-00599C?style=for-the-badge">

<img src="https://img.shields.io/badge/Research%20Status-Completed-2E7D32?style=for-the-badge">

<img src="https://img.shields.io/badge/Publication-2026-0A7EA4?style=for-the-badge">

<img src="https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge">

</p>

<br>

### 🐘 Human TP53 ↔ Elephant TP53

**Comparative genomics • Cancer biology • Evolutionary bioinformatics**

</div>

---

# 📑 Research Navigation

- [🔎 Research at a Glance](#-research-at-a-glance)
- [📄 Abstract](#-abstract)
- [🧠 Biological Rationale](#-biological-rationale)
- [🎯 Aim & Objectives](#-aim--objectives)
- [🗂️ Data Sources](#️-data-sources)
- [🔬 Computational Methodology](#-computational-methodology)
- [🧰 Computational Tools](#-computational-tools)
- [📊 Key Findings](#-key-findings)
- [🐘 EleProtect](#-eleprotect)
- [📷 Research Visualisation](#-research-visualisation)
- [📂 Repository Structure](#-repository-structure)
- [♻️ Reproducibility](#️-reproducibility)
- [📚 Publication](#-publication)
- [🎓 Academic Research Context](#-academic-research-context)
- [🔭 Future Research](#-future-research)
- [⚠️ Scope & Responsible Interpretation](#️-scope--responsible-interpretation)
- [🙏 Acknowledgements](#-acknowledgements)
- [📖 Citation](#-citation)
- [👩‍🔬 Author](#-author)

---

# 🔎 Research at a Glance

| Category | Details |
|----------|---------|
| 🧬 **Research Area** | Comparative Genomics • Cancer Biology • Evolutionary Bioinformatics |
| 🧬 **Primary Gene** | **TP53 / p53** |
| 🐘 **Comparative System** | Human • African Elephant • Asian Elephant |
| 🎯 **Focal Hotspots** | R175 • G245 • R248 • R249 • R273 • R282 |
| 🔬 **Approach** | Protein-level multiple sequence alignment • hotspot mapping • conservation analysis |
| 🧪 **Core Tools** | NCBI • UniProt • MAFFT • Jalview • MEGA • Python • Google Colab |
| 📊 **Research Status** | ✅ Completed |
| 📚 **Publication Status** | ✅ Published |
| 🐘 **Companion Application** | EleProtect |
| 🎓 **Academic Context** | Postgraduate Bioinformatics Research |

---

# 📄 Abstract

This research presents a comparative in-silico analysis of **TP53 mutation hotspots** in humans and elephants, with particular emphasis on sequence conservation and variation across canonical and retrogene-derived elephant TP53 sequences.

The study was motivated by the evolutionary question underlying **Peto's paradox**: despite their large body size and long lifespan, elephants do not exhibit cancer incidence proportional to the number of cells or years at risk.

Major human TP53 mutation hotspots were mapped against elephant TP53 protein sequences to investigate whether residues frequently altered in human cancer show conservation across elephant lineages and whether corresponding retrogene sequences display distinct patterns of variation.

The analysis integrates publicly available biological sequence resources with established computational biology workflows, including multiple sequence alignment, hotspot mapping, conservation assessment, and comparative sequence analysis.

The project demonstrates how transparent comparative genomics can be used to investigate evolutionary constraints surrounding tumour-suppressor genes and generate hypotheses for future experimental research.

---

# 🧠 Biological Rationale

**TP53** encodes the tumour-suppressor protein p53, a central regulator of genomic stability, DNA-damage response, cell-cycle control, senescence, and apoptosis.

Recurrent mutations in TP53 occur at functionally important positions in human cancers. These recurrent hotspots provide a useful framework for investigating whether evolutionarily conserved residues remain constrained across species.

Elephants represent a particularly interesting comparative system because of their large body size and relatively long lifespan.

The relationship between body size, lifespan, cell number, and cancer incidence is not proportional in the way that simple scaling arguments would predict. This evolutionary observation is commonly discussed in the context of **Peto's paradox**.

Rather than assuming that sequence conservation directly explains cancer resistance, this study uses comparative computational analysis to identify patterns of conservation and divergence that may inform future functional and experimental investigation.

---

# 🎯 Aim & Objectives

## Aim

To investigate the evolutionary conservation and sequence-level variation of major human TP53 mutation hotspots across elephant TP53 sequences using comparative in-silico bioinformatics approaches.

## Objectives

- 🧬 Identify recurrent human TP53 cancer mutation hotspots.
- 🗂️ Retrieve human and elephant TP53 sequences from publicly available biological databases.
- 🔬 Perform protein-level multiple sequence alignment.
- 🎯 Map human TP53 mutation hotspots onto elephant TP53 sequences.
- 📊 Evaluate conservation and amino-acid variation at corresponding positions.
- 🐘 Compare canonical elephant TP53 sequences with TP53 retrogene-derived sequences.
- 🧪 Explore sequence-level patterns relevant to evolutionary cancer biology.
- ♻️ Establish a transparent computational framework that can be inspected and extended.

---

# 🗂️ Data Sources

The study uses publicly accessible biological sequence and annotation resources.

| Resource | Purpose |
|----------|---------|
| 🧬 **NCBI** | Sequence retrieval and biological database resources |
| 🧬 **UniProt** | Protein sequence and annotation resources |
| 🧬 **Human TP53 mutation data** | Identification of recurrent mutation hotspots |
| 🐘 **Elephant TP53 sequences** | Comparative sequence analysis |
| 🧬 **Elephant TP53 retrogenes** | Assessment of sequence-level variation |

### Data Principle

The project is based on **publicly accessible biological data** and does not require patient-identifiable genomic information.

---

# 🔬 Computational Methodology

The overall computational workflow can be represented as:

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
 Multiple Sequence Alignment
            │
            ▼
     Hotspot Position Mapping
            │
            ▼
 Conservation & Variation Analysis
            │
            ▼
 Canonical vs Retrogene Comparison
            │
            ▼
  Descriptive Hotspot Prioritisation
            │
            ▼
     Biological Interpretation
