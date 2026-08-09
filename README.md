# 🧬 Elephant TP53 Hotspot Mapping

## 🧑 Human TP53 ↔ 🐘 Elephant TP53

### Comparative In-Silico Analysis of TP53 Mutation Hotspots Between Humans and Elephants

**Comparative Genomics · Cancer Biology · Evolutionary Bioinformatics · Computational Biology · Research Software · Reproducible Research**

<p align="center">

<a href="https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping">
<img src="https://img.shields.io/badge/Research-Comparative%20Genomics-1f6feb" alt="Research">
</a>
<a href="https://www.python.org/">
<img src="https://img.shields.io/badge/Python-3.x-3776AB" alt="Python">
</a>
<a href="https://streamlit.io/">
<img src="https://img.shields.io/badge/EleProtect-Streamlit-FF4B4B" alt="EleProtect">
</a>
<a href="LICENSE">
<img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License">
</a>
<a href="https://eccb2026.org/">
<img src="https://img.shields.io/badge/ECCB%202026-Poster%20Presentation-orange" alt="ECCB 2026">
</a>

</p>

---

# 📑 Research Navigation

- [🔎 Research at a Glance](#-research-at-a-glance)
- [📄 Abstract](#-abstract)
- [🧠 Biological Rationale](#-biological-rationale)
- [🧑 Human TP53 Reference](#-human-tp53-reference)
- [🎯 TP53 Mutation-Associated Hotspots](#-tp53-mutation-associated-hotspots)
- [🐘 Why Elephants?](#-why-elephants)
- [🧬 Research Framework](#-research-framework)
- [🎯 Aim & Objectives](#-aim--objectives)
- [🔬 Computational Methodology](#-computational-methodology)
- [🗂️ Data Sources](#️-data-sources)
- [🐘 Elephant Genomic Reference](#-elephant-genomic-reference)
- [📁 Repository Structure](#-repository-structure)
- [🧭 Recommended Path for a Researcher / PI](#-recommended-path-for-a-researcher--pi)
- [📊 Research Outputs](#-research-outputs)
- [♻️ Reproducibility](#️-reproducibility)
- [🐍 Reusable Computational Code](#-reusable-computational-code)
- [📓 Main Research Notebooks](#-main-research-notebooks)
- [💻 EleProtect v2.0](#-eleprotect-v20)
- [⚠️ Scientific Scope & Limitations](#️-scientific-scope--limitations)
- [📚 Published Research](#-published-research)
- [📄 Research Square](#-research-square)
- [🌍 ECCB 2026](#-eccb-2026)
- [🎓 Academic Context](#-academic-context)
- [🙏 Acknowledgements](#-acknowledgements)
- [👥 Authors & Contributors](#-authors--contributors)
- [🏆 Researcher Profile](#-researcher-profile)
- [🔗 Connect With the Researcher](#-connect-with-the-researcher)
- [📖 Citation](#-citation)
- [📜 License](#-license)
- [🌐 Research Perspective](#-research-perspective)
- [🧭 Research Philosophy](#-research-philosophy)
- [🚀 Future Research Direction](#-future-research-direction)
- [🔗 Quick Access](#-quick-access)

---

# 🔎 Research at a Glance

| Category | Details |
|---|---|
| 🧬 **Research Area** | Comparative Genomics · Cancer Biology · Evolutionary Bioinformatics |
| 🧬 **Primary Gene** | **TP53 / p53** |
| 🧑 **Human Reference System** | *Homo sapiens* |
| 🐘 **Comparative System** | *Elephas maximus* and TP53-related elephant sequences |
| 🎯 **Focal Hotspots** | R175 · G245 · R248 · R249 · R273 · R282 |
| 🔬 **Approach** | Protein-level multiple sequence alignment · hotspot mapping · conservation analysis · phylogenetic analysis |
| 🧪 **Core Resources** | NCBI · UniProt · NCBI BLAST · MAFFT · Jalview · MEGA |
| 💻 **Computational Environment** | Python · Google Colab |
| 🤖 **Exploratory Analysis** | Sequence-derived feature analysis · computational prioritization · exploratory machine learning |
| 💻 **Research Software** | EleProtect v2.0 |
| 📊 **Research Status** | Completed |
| 📚 **Publication Status** | Published |
| 🌍 **Conference Achievement** | ECCB 2026 Poster Presentation |
| 🎓 **Academic Context** | Postgraduate Bioinformatics Research |

---

# 📄 Abstract

This research presents a comparative in-silico analysis of **TP53 mutation-associated hotspots** in humans and elephants, with particular emphasis on sequence conservation and variation across canonical and TP53-related elephant sequences.

The study was motivated by the evolutionary question underlying **Peto's paradox**: despite their large body size and long lifespan, elephants do not exhibit cancer incidence proportional to the number of cells or years at risk.

Major human TP53 mutation-associated hotspots were mapped against elephant TP53-related protein sequences to investigate whether residues frequently altered in human cancer show conservation across elephant sequences and whether corresponding TP53-related sequences display distinct patterns of variation.

The analysis integrates publicly available biological sequence resources with established computational biology workflows, including similarity-based sequence identification, sequence curation, multiple sequence alignment, hotspot mapping, conservation assessment, phylogenetic analysis, sequence-derived feature analysis, and exploratory computational prioritization.

The project demonstrates how transparent comparative genomics can be used to investigate evolutionary constraints surrounding tumour-suppressor genes and generate hypotheses for future structural, functional, and experimental research.

> **Scientific interpretation:** The computational observations generated by this project are intended to identify sequence-level patterns and research hypotheses. They do not independently establish a complete molecular explanation for elephant cancer resistance or Peto's paradox.

---

# 🧠 Biological Rationale

**TP53** encodes the tumour-suppressor protein p53, a central regulator of genomic stability, DNA-damage response, cell-cycle control, senescence, and apoptosis.

Recurrent mutations in TP53 occur at functionally important positions in human cancers. These mutation-associated hotspots provide a useful framework for investigating whether evolutionarily constrained residues remain conserved across species.

The comparative analysis of human and elephant TP53-related sequences provides an evolutionary framework for examining sequence conservation and divergence at positions associated with human TP53 cancer mutations.

The relationship between body size, lifespan, cell number, and cancer incidence is commonly discussed in the context of **Peto's paradox**.

Rather than assuming that sequence conservation directly explains cancer resistance, this study uses comparative computational analysis to identify sequence-level patterns that may inform future biological and experimental investigation.

---

# 🧑 Human TP53 Reference

The human TP53 protein provides the primary reference framework for the comparative analysis.

| Attribute | Details |
|---|---|
| **Species** | *Homo sapiens* |
| **Gene** | **TP53** |
| **Protein** | p53 |
| **UniProt Accession** | **P04637** |
| **Sequence Type** | Canonical human TP53 protein |
| **Role in Analysis** | Reference sequence for hotspot mapping and comparative sequence analysis |

🔗 **[UniProt P04637 — Human TP53](https://www.uniprot.org/uniprotkb/P04637/entry)**

---

# 🎯 TP53 Mutation-Associated Hotspots

The comparative framework focuses on recurrent human TP53 mutation-associated positions:

| Hotspot | Human TP53 Position |
|---|---:|
| **R175** | 175 |
| **G245** | 245 |
| **R248** | 248 |
| **R249** | 249 |
| **R273** | 273 |
| **R282** | 282 |

These positions provide defined human reference coordinates for cross-species residue-level comparison following sequence alignment.

The workflow does not assume that conservation at these positions is itself evidence of cancer resistance. Conservation and amino-acid variation are treated as computational observations that can support subsequent biological investigation.

---

# 🐘 Why Elephants?

Elephants provide an important comparative system in evolutionary cancer biology because they combine:

```text
Large body size
       +
Long lifespan
       +
Large number of cells
       ↓
A major evolutionary cancer-biology question
```

The relationship between body size, lifespan, cell number, and cancer incidence is commonly discussed in the context of **Peto's paradox**.

This project therefore uses elephants as a comparative evolutionary system to investigate TP53 sequence conservation and variation.

> **Important:** The research does not claim that elephant cancer resistance is explained by TP53 sequence conservation alone.

Instead, the computational objective is to identify sequence-level patterns that may inform future biological investigation.

---

# 🧬 Research Framework

The overall analytical logic is:

```text
🧑 Human TP53
       ↓
🎯 Human Cancer-Associated Hotspots
       ↓
🐘 Elephant Protein Resources
       ↓
🔎 Similarity-Based Sequence Identification
       ↓
🧹 Sequence Curation
       ↓
🧬 Multiple Sequence Alignment
       ↓
🎯 Hotspot Position Mapping
       ↓
📊 Conservation & Amino-Acid Variation
       ↓
🐘 Canonical vs TP53-Related Comparison
       ↓
🌳 Phylogenetic Analysis
       ↓
📊 Sequence-Derived Feature Analysis
       ↓
🤖 Exploratory Computational Prioritization
       ↓
🧠 Biological Interpretation
       ↓
🔬 Hypothesis Generation
```

---

# 🎯 Aim & Objectives

## Aim

To investigate the evolutionary conservation and sequence-level variation of major human TP53 mutation-associated hotspots across elephant TP53-related sequences using comparative in-silico bioinformatics approaches.

## Objectives

- Identify recurrent human TP53 mutation-associated hotspots.
- Retrieve human and elephant TP53-related sequence resources.
- Identify TP53-related elephant sequences using similarity-based analysis.
- Curate sequences for comparative analysis.
- Perform protein-level multiple sequence alignment.
- Map human TP53 hotspots onto elephant sequences.
- Evaluate conservation and amino-acid variation.
- Compare canonical elephant TP53 with TP53-related sequences.
- Examine evolutionary relationships among analysed sequences.
- Develop exploratory sequence-derived computational prioritization.
- Investigate exploratory machine-learning approaches.
- Preserve the computational workflow as a transparent and reproducible research resource.
- Develop **EleProtect v2.0**, an interactive research application associated with the analysis.

---

# 🔬 Computational Methodology

The project follows a multi-stage comparative bioinformatics workflow.

## 1. Sequence Identification

Human TP53 is used as a reference for identifying TP53-related sequences within elephant protein resources.

Similarity-based analysis is used to identify candidate sequences for subsequent evaluation.

## 2. Sequence Curation

Candidate sequences are inspected and organized into appropriate datasets for comparative analysis.

The repository distinguishes between:

- Raw biological resources
- Processed sequences
- Alignment-ready datasets
- Analytical outputs

## 3. Multiple Sequence Alignment

Protein sequences are aligned to establish residue-level correspondence between human TP53 and elephant TP53-related sequences.

Primary computational resources include:

- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [Jalview](https://www.jalview.org/)
- [MEGA](https://www.megasoftware.net/)

## 4. Hotspot Mapping

Human TP53 mutation-associated positions are mapped onto aligned sequences.

```text
Human TP53 hotspot
        ↓
Reference residue position
        ↓
Alignment-aware position
        ↓
Elephant corresponding residue
        ↓
Conservation / substitution
        ↓
Comparative interpretation
```

## 5. Conservation Analysis

The analysis evaluates sequence conservation and amino-acid variation at positions corresponding to selected human TP53 hotspots.

## 6. Canonical vs TP53-Related Sequence Comparison

Canonical elephant TP53 and TP53-related / retrogene-derived sequences are compared to investigate sequence-level variation.

## 7. Phylogenetic Analysis

Phylogenetic analysis provides evolutionary context for the analysed sequences.

## 8. Sequence-Derived Feature Analysis

Sequence-derived characteristics are calculated for analysed sequences, including composition-related features and similarity measures.

## 9. Exploratory Computational Prioritization

Sequence-derived features are used in exploratory computational and machine-learning workflows.

> **Important:** These analyses are exploratory and are not presented as clinically validated prediction models.

---

# 🗂️ Data Sources

The project uses publicly accessible biological sequence resources.

| Resource | Purpose |
|---|---|
| [NCBI](https://www.ncbi.nlm.nih.gov/) | Biological sequence and genomic resources |
| [UniProt](https://www.uniprot.org/) | Protein sequence and annotation resources |
| [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) | Similarity-based sequence identification |
| Elephant genomic resources | Comparative elephant sequence analysis |
| Human TP53 resources | Reference sequence and hotspot analysis |

---

# 🐘 Elephant Genomic Reference

| Attribute | Details |
|---|---|
| **Species** | *Elephas maximus* |
| **Assembly** | `GCA_024166365.1` |
| **Assembly Name** | mEleMax1 primary haplotype |
| **Source** | NCBI Assembly |

🔗 **[NCBI Assembly — GCA_024166365.1](https://www.ncbi.nlm.nih.gov/assembly/)**

For complete sequence provenance, accession information, and source documentation:

👉 **[Read `docs/provenance.md`](docs/provenance.md)**

---

# 📁 Repository Structure

```text
Elephant-TP53-Hotspot-Mapping/
│
├── data/
│   ├── raw/
│   │   ├── Blast/
│   │   ├── ncbi_dataset/
│   │   ├── human_tp53.fasta
│   │   ├── elephant_tp53_canonical.fasta
│   │   ├── elephant_tp53_RTG9.fasta
│   │   ├── elephant_tp53_RTG9_translate.fasta
│   │   ├── african_tp53_hits.fasta
│   │   ├── asian_tp53_hits.fasta
│   │   └── provenance.txt
│   │
│   ├── processed/
│   │   ├── TP53_clean.fasta
│   │   ├── TP53_all_sequences.fasta
│   │   ├── human_elephant_tp53_pair.fasta
│   │   └── human_elephant_tp53_retrogene_comparison.fasta
│   │
│   └── Database/
│       ├── README.md
│       ├── tp53_elephant_database.csv
│       └── tp53_elephant_database.json
│
├── notebooks/
│   ├── TP53_Comparative_Analysis.ipynb
│   ├── EleProtect_Feature_Prioritization.ipynb
│   └── TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb
│
├── scripts/
│   └── TP53_Comparative_Analysis.py
│
├── results/
│   ├── MEGA/
│   ├── ML/
│   ├── MSA/
│   ├── MSA_1and2/
│   ├── MSA_3/
│   ├── phylogeny/
│   └── README.md
│
├── figures/
│
├── docs/
│   ├── methodology.md
│   ├── provenance.md
│   ├── reproducibility.md
│   ├── interpretation_and_limitations.md
│   └── Research_Proposal.pdf
│
├── EleProtect_App/
│
├── CITATION.cff
├── LICENSE
├── README.md
└── .gitignore
```

---

# 🧭 Recommended Path for a Researcher / PI

If you are reviewing this repository for the first time:

### 01 — 🔬 Understand the Research

Read this README to understand the scientific question, biological rationale, computational framework, and research outputs.

👉 **[Start with the main README](README.md)**

### 02 — 🧪 Inspect the Methodology

Review the detailed computational methodology.

👉 **[docs/methodology.md](docs/methodology.md)**

### 03 — 🧾 Verify Data Provenance

Inspect biological sources, accession information, and sequence provenance.

👉 **[docs/provenance.md](docs/provenance.md)**

### 04 — ♻️ Inspect Reproducibility

Review the reproducibility framework.

👉 **[docs/reproducibility.md](docs/reproducibility.md)**

### 05 — 📓 Open the Main Analysis

Explore the primary comparative analysis.

👉 **[TP53_Comparative_Analysis.ipynb](notebooks/TP53_Comparative_Analysis.ipynb)**

### 06 — 📊 Inspect Computational Outputs

Explore the generated computational outputs.

👉 **[results/](results/)**

### 07 — 🖼️ Examine Figures

Explore research visualizations.

👉 **[figures/](figures/)**

### 08 — 🤖 Explore Feature Prioritization

Review the exploratory sequence-feature analysis.

👉 **[EleProtect_Feature_Prioritization.ipynb](notebooks/EleProtect_Feature_Prioritization.ipynb)**

### 09 — 💻 Explore the Research Application

Explore the associated research software.

👉 **[EleProtect_App/](EleProtect_App/)**

### 10 — ⚠️ Review Scientific Limitations

Review the interpretation boundaries of the computational findings.

👉 **[docs/interpretation_and_limitations.md](docs/interpretation_and_limitations.md)**

### 11 — 📚 Explore Associated Research

See the publication, Research Square work, and ECCB research direction below.

---

# 📊 Research Outputs

The repository preserves computational outputs generated during the study.

## 🧬 Multiple Sequence Alignment

- [results/MSA/](results/MSA/)
- [results/MSA_1and2/](results/MSA_1and2/)
- [results/MSA_3/](results/MSA_3/)

## 🌳 Phylogenetic Analysis

- [results/phylogeny/](results/phylogeny/)
- [results/MEGA/](results/MEGA/)

## 🤖 Computational / ML Outputs

- [results/ML/](results/ML/)

## 🖼️ Research Figures

- [figures/](figures/)

## 📋 Results Documentation

- [results/README.md](results/README.md)

---

# ♻️ Reproducibility

A central objective of this repository is to maintain a traceable relationship between the biological source data, computational analysis, and resulting observations.

```text
🧬 Biological Source
        ↓
🧾 Data Provenance
        ↓
📂 Raw Data
        ↓
⚙️ Processing
        ↓
🔬 Analysis
        ↓
📊 Results
        ↓
🖼️ Figures
        ↓
🧠 Interpretation
```

Detailed documentation:

- 🧪 **[Methodology](docs/methodology.md)**
- 🧾 **[Data Provenance](docs/provenance.md)**
- ♻️ **[Reproducibility](docs/reproducibility.md)**
- ⚠️ **[Interpretation & Limitations](docs/interpretation_and_limitations.md)**

---

# 🐍 Reusable Computational Code

Reusable computational implementations are maintained under:

👉 **[scripts/](scripts/)**

The principal comparative implementation is:

👉 **[scripts/TP53_Comparative_Analysis.py](scripts/TP53_Comparative_Analysis.py)**

The repository separates reusable computational implementation from interactive analysis:

```text
scripts/
    ↓
Reusable computational implementation

notebooks/
    ↓
Interactive analysis, exploration and interpretation
```

---

# 📓 Main Research Notebooks

## 1. 🧬 TP53 Comparative Analysis

Primary computational workflow for the comparative TP53 analysis.

👉 **[Open TP53_Comparative_Analysis.ipynb](notebooks/TP53_Comparative_Analysis.ipynb)**

## 2. 🤖 EleProtect Feature Prioritization

Feature extraction and exploratory computational prioritization.

👉 **[Open EleProtect_Feature_Prioritization.ipynb](notebooks/EleProtect_Feature_Prioritization.ipynb)**

## 3. 🧠 Predictive Modelling

Exploratory computational and predictive-modelling workflow.

👉 **[Open TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb](notebooks/TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb)**

---

# 💻 EleProtect v2.0

## Interactive Research Application

**EleProtect v2.0** is the software component associated with this research project.

It provides a Streamlit-based interface for selected sequence-analysis and exploratory computational-prioritization workflows.

### Features

- 🧬 DNA / protein sequence input
- 🔎 Sequence processing
- 🎯 TP53-oriented analysis
- 📊 Feature extraction
- 🤖 Exploratory ML scoring
- 📁 CSV export
- 🌐 Interactive Streamlit interface

### Explore EleProtect

👉 **[Open EleProtect_App](EleProtect_App/)**

👉 **[Read EleProtect_App/README.md](EleProtect_App/README.md)**

### Local Execution

```bash
cd EleProtect_App
pip install -r requirements.txt
streamlit run app.py
```

> ⚠️ **Research-use statement:** EleProtect is a research prototype and is not intended to function as a clinical diagnostic, prognostic, or therapeutic decision-support system.

---

# ⚠️ Scientific Scope & Limitations

This repository contains computational research and should be interpreted within the scope of the available sequence-level evidence.

### The analysis can provide evidence concerning:

- Sequence similarity
- Sequence conservation
- Amino-acid variation
- TP53-related sequence relationships
- Phylogenetic context
- Sequence-derived computational prioritization

### The analysis does not independently establish:

- The complete molecular basis of elephant cancer resistance
- A complete mechanistic explanation of Peto's paradox
- Functional activity of individual TP53-related retrogenes
- Clinical cancer risk
- Therapeutic response
- Diagnostic utility
- Experimentally validated biological mechanisms

The intended scientific progression is:

```text
Computational Observation
        ↓
Evolutionary Interpretation
        ↓
Hypothesis Generation
        ↓
Structural / Functional Investigation
        ↓
Experimental Validation
```

👉 **[Read Interpretation & Limitations](docs/interpretation_and_limitations.md)**

---

# 📚 Published Research

## Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants

**Rawat, R. R., Nadar, S., & Uppal, G. K. (2026).**

*Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants: A Responsible Bioinformatics Innovation Contributing to Cancer Research.*

*Ideathon on Vikshit Bharat: Ideas, Innovation and Impact.*

**Chapter 11, pp. 123–138.**

### 🔗 Publication

👉 **[DOI: 10.25215/9141002199](https://doi.org/10.25215/9141002199)**

### 🧾 Machine-Readable Citation

👉 **[CITATION.cff](CITATION.cff)**

If you use this repository, its methodology, computational framework, figures, research software, or associated research framework, please cite the associated publication.

---

# 📄 Research Square

The broader TP53 research direction extends this comparative work toward investigation of **TP53 hotspot conservation and functional constraint across mammalian species**.

### 🔗 Research Square Publication

👉 **[DOI: 10.21203/rs.3.rs-9299199/v1](https://doi.org/10.21203/rs.3.rs-9299199/v1)**

---

# 🌍 ECCB 2026

## European Conference on Computational Biology

The broader TP53 research programme has been **accepted for poster presentation at ECCB 2026**.

### Poster Title

**Evolutionary Conservation and Functional Constraint of TP53 Mutation Hotspots Across Mammalian Species**

### Conference

**ECCB 2026 — European Conference on Computational Biology**

📍 **Geneva, Switzerland**

📅 **31 August – 4 September 2026**

🎓 **Poster Presentation**

🔗 **[ECCB 2026 — European Conference on Computational Biology](https://eccb2026.org/)**

### 🔗 How This Repository Connects to ECCB

The ECCB research represents a broader continuation of the TP53 research direction developed through this project.

```text
🧑 Human TP53 Mutation Hotspots
              ↓
🐘 Human–Elephant Comparative Analysis
              ↓
🧬 TP53 Conservation & Sequence Variation
              ↓
🔬 Canonical / TP53-Related Sequence Comparison
              ↓
📚 Published Research
              ↓
🌳 Broader Mammalian TP53 Analysis
              ↓
📄 Research Square
              ↓
🌍 ECCB 2026 Poster Presentation
```

This repository therefore represents one component of a continuing research trajectory in:

**Comparative Genomics · Evolutionary Cancer Biology · Computational Biology**

---

# 🎓 Academic Context

This research was developed within the postgraduate academic and research context of:

### Department of Bioinformatics

**Guru Nanak Khalsa College of Arts, Commerce and Science**

**University of Mumbai**

The project represents a computational research direction connecting:

```text
Bioinformatics
      ↓
Comparative Genomics
      ↓
Cancer Biology
      ↓
Evolutionary Analysis
      ↓
Computational Research
      ↓
Research Software
```

---

# 🙏 Acknowledgements

The project acknowledges the academic environment associated with:

**Department of Bioinformatics**  
**Guru Nanak Khalsa College of Arts, Commerce and Science**  
**University of Mumbai**

The project also acknowledges the publicly available biological resources and computational tools that supported the analysis:

- [NCBI](https://www.ncbi.nlm.nih.gov/)
- [UniProt](https://www.uniprot.org/)
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [Jalview](https://www.jalview.org/)
- [MEGA](https://www.megasoftware.net/)
- [Python](https://www.python.org/)
- [Google Colab](https://colab.research.google.com/)
- [scikit-learn](https://scikit-learn.org/)
- [Streamlit](https://streamlit.io/)

---

# 👥 Authors & Contributors

## 👩‍🔬 Ritika Rajendra Rawat

**Primary Researcher / Author**

Contributed to the computational research framework, comparative analysis, sequence analysis, research development, and associated research outputs.

## 🧑‍🔬 Sermarani Nadar

**Co-author**

Co-author of the associated published research on comparative TP53 hotspot mapping in elephants.

## 🧑‍🔬 Gursimran Kaur Uppal

**Co-author**

Co-author of the associated published research on comparative TP53 hotspot mapping in elephants.

---

# 🏆 Researcher Profile

## 👩‍🔬 Ritika Rajendra Rawat

**MSc Bioinformatics**

Computational biology researcher working at the intersection of:

- 🧬 Computational Genomics
- 🧬 Comparative Genomics
- 🧪 Cancer Bioinformatics
- 🌳 Evolutionary Bioinformatics
- 🧬 Protein & Sequence Analysis
- 🤖 Machine Learning for Biological Data
- 💻 Research Software
- ♻️ Reproducible Computational Research

### 🔬 Research Direction

My research interests focus on applying computational approaches to biological and evolutionary questions, particularly:

> **How can sequence-level computational evidence reveal evolutionary constraints relevant to cancer biology?**

This project represents an effort to build an **open, inspectable, and reproducible computational research framework** rather than treating computational results as isolated findings.

---

# 🏅 Research Highlights & Achievements

| Achievement | Details |
|---|---|
| 🎓 **Academic** | MSc Bioinformatics |
| 🧬 **Research Focus** | Comparative TP53 and evolutionary cancer biology |
| 📚 **Publication** | Co-authored published TP53 comparative research |
| 📄 **Research Square** | Broader mammalian TP53 research |
| 🌍 **ECCB 2026** | Accepted poster presentation |
| 💻 **Research Software** | Developer of EleProtect v2.0 |
| ♻️ **Open Research** | Reproducible GitHub research repository |
| 🧬 **Research Output** | Comparative TP53 hotspot analysis |

---

# 🔗 Connect With the Researcher

## 💻 GitHub

**Rita1791**

👉 **[GitHub Profile](https://github.com/Rita1791)**

👉 **[Elephant TP53 Hotspot Mapping Repository](https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping)**

---

## 💼 LinkedIn

**Ritika Rajendra Rawat**

👉 **[Connect with Ritika on LinkedIn](https://in.linkedin.com/in/ritika-rawat-551107219)**

Academic and professional networking:

- Research collaboration
- Computational biology discussions
- Scientific networking
- PhD-related academic communication

---

## 🔬 ResearchGate

**Ritika Rawat**

👉 **[View Ritika's ResearchGate Profile](https://www.researchgate.net/profile/Ritika-Rawat-10)**

ResearchGate provides access to the researcher's broader research and publication profile.

---

## 📧 Research / Academic Email

**ritika.rawat27@outlook.com**

👉 **[Contact Ritika by Email](mailto:ritika.rawat27@outlook.com)**

Suitable for:

- Academic collaboration
- Research discussions
- Computational biology projects
- PhD-related academic communication
- Scientific enquiries

---

# 📖 Citation

If you use this repository, its methodology, computational framework, figures, research software, or associated research framework, please cite the associated publication.

### Associated Publication

> **Rawat, R. R., Nadar, S., & Uppal, G. K. (2026).**  
> *Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants: A Responsible Bioinformatics Innovation Contributing to Cancer Research.*  
> *Ideathon on Vikshit Bharat: Ideas, Innovation and Impact.*  
> Chapter 11, pp. 123–138.  
> DOI: [10.25215/9141002199](https://doi.org/10.25215/9141002199)

### Machine-Readable Citation

👉 **[CITATION.cff](CITATION.cff)**

---

# 📜 License

This repository is released under the **MIT License**.

👉 **[View LICENSE](LICENSE)**

---

# 🌐 Research Perspective

This project brings together:

```text
        🧬 Comparative Genomics
                  +
        🧪 Cancer Biology
                  +
        🌳 Evolutionary Biology
                  +
        💻 Bioinformatics
                  +
        📊 Sequence Analysis
                  +
        🤖 Computational Modelling
                  +
        🧰 Research Software
                  +
        ♻️ Reproducible Research
```

The central principle of this repository is:

> **Use computational evidence to identify evolutionary patterns, preserve the analytical path behind those observations, and translate computational findings into hypotheses that can ultimately be investigated through structural, functional, and experimental research.**

---

# 🧭 Research Philosophy

The project is designed around a simple principle:

```text
Transparent Data
      +
Traceable Analysis
      +
Reproducible Computation
      +
Responsible Interpretation
      =
Open Computational Research
```

The repository therefore preserves not only selected results, but also the computational pathway through which those results were generated.

---

# 🚀 Future Research Direction

Potential future extensions include:

- Structural analysis of conserved and divergent TP53 regions
- Comparative analysis across additional mammalian species
- Functional investigation of TP53-related sequences
- Integration of structural and evolutionary features
- Expanded comparative cancer-genomics datasets
- More rigorous benchmarking of computational prioritization methods
- Experimental validation of biologically relevant hypotheses

The intended progression is:

```text
Computational Analysis
        ↓
Comparative Observation
        ↓
Evolutionary Hypothesis
        ↓
Structural Investigation
        ↓
Functional Investigation
        ↓
Experimental Validation
```

---

# 🔗 Quick Access

| Research Resource | Link |
|---|---|
| 🧬 **GitHub Repository** | [Elephant-TP53-Hotspot-Mapping](https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping) |
| 🧪 **Methodology** | [docs/methodology.md](docs/methodology.md) |
| 🧾 **Provenance** | [docs/provenance.md](docs/provenance.md) |
| ♻️ **Reproducibility** | [docs/reproducibility.md](docs/reproducibility.md) |
| ⚠️ **Interpretation & Limitations** | [docs/interpretation_and_limitations.md](docs/interpretation_and_limitations.md) |
| 📊 **Results** | [results/](results/) |
| 🖼️ **Figures** | [figures/](figures/) |
| 📓 **Main Analysis** | [TP53_Comparative_Analysis.ipynb](notebooks/TP53_Comparative_Analysis.ipynb) |
| 🤖 **Feature Prioritization** | [EleProtect_Feature_Prioritization.ipynb](notebooks/EleProtect_Feature_Prioritization.ipynb) |
| 💻 **EleProtect** | [EleProtect_App/](EleProtect_App/) |
| 📚 **Published Research** | [DOI: 10.25215/9141002199](https://doi.org/10.25215/9141002199) |
| 📄 **Research Square** | [DOI: 10.21203/rs.3.rs-9299199/v1](https://doi.org/10.21203/rs.3.rs-9299199/v1) |
| 🌍 **ECCB 2026** | [European Conference on Computational Biology](https://eccb2026.org/) |
| 🔬 **ResearchGate** | [Ritika Rawat](https://www.researchgate.net/profile/Ritika-Rawat-10) |
| 💼 **LinkedIn** | [Ritika Rajendra Rawat](https://in.linkedin.com/in/ritika-rawat-551107219) |
| 📧 **Email** | [ritika.rawat27@outlook.com](mailto:ritika.rawat27@outlook.com) |
| 📜 **License** | [MIT License](LICENSE) |
| 🧾 **Citation File** | [CITATION.cff](CITATION.cff) |

---

<p align="center">

### 🧑 Human TP53 ↔ 🐘 Elephant TP53

**Comparative Genomics · Evolutionary Cancer Biology · Computational Biology**

<br>

<a href="https://doi.org/10.25215/9141002199">📚 Published Research</a>
&nbsp; · &nbsp;
<a href="https://doi.org/10.21203/rs.3.rs-9299199/v1">📄 Research Square</a>
&nbsp; · &nbsp;
<a href="https://eccb2026.org/">🌍 ECCB 2026</a>
&nbsp; · &nbsp;
<a href="EleProtect_App/">💻 EleProtect v2.0</a>

<br><br>

### 👩‍🔬 Ritika Rajendra Rawat

**MSc Bioinformatics**

<br>

<a href="https://github.com/Rita1791">GitHub</a>
&nbsp; · &nbsp;
<a href="https://in.linkedin.com/in/ritika-rawat-551107219">LinkedIn</a>
&nbsp; · &nbsp;
<a href="https://www.researchgate.net/profile/Ritika-Rawat-10">ResearchGate</a>
&nbsp; · &nbsp;
<a href="mailto:ritika.rawat27@outlook.com">Email</a>

</p>
