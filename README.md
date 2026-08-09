Yes. Below is the **complete root `README.md` as one single Markdown file**. Copy everything inside the code block directly into the repository's **root `README.md`**.

````markdown
# 🧬 Elephant TP53 Hotspot Mapping

## Comparative In-Silico Analysis of TP53 Mutation Hotspots Between Humans and Elephants

<p align="center">

**Comparative Genomics · Cancer Biology · Evolutionary Bioinformatics · Computational Biology**

<br><br>

🧬 **Human TP53** &nbsp; ↔ &nbsp; 🐘 **Elephant TP53**

</p>

<p align="center">

[![Research](https://img.shields.io/badge/Research-Comparative%20Genomics-1f6feb)](https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping)
[![Python](https://img.shields.io/badge/Python-3.x-3776AB)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/EleProtect-Streamlit-FF4B4B)](https://streamlit.io/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![ECCB 2026](https://img.shields.io/badge/ECCB%202026-Poster%20Presentation-orange)](https://eccb2026.org/)

</p>

---

# 🔎 Research at a Glance

**Elephant TP53 Hotspot Mapping** is a comparative bioinformatics research
project investigating the **evolutionary conservation and sequence-level
variation of major human TP53 mutation-associated hotspots across elephant
TP53-related sequences**.

The project was developed in the context of **evolutionary cancer biology**
and the biological question commonly associated with **Peto's paradox**:
cancer incidence does not simply increase proportionally with body size,
number of cells, and lifespan across species.

Elephants provide an important comparative system for investigating this
question because of their large body size and long lifespan. This project
focuses on **TP53**, a major tumour-suppressor gene involved in genomic
stability, DNA-damage response, cell-cycle regulation, apoptosis, and cellular
stress responses.

Rather than assuming that one gene or one sequence feature explains elephant
cancer resistance, this research uses comparative sequence analysis to
investigate whether **human TP53 mutation-associated hotspot residues are
conserved, substituted, or otherwise variable across elephant TP53-related
sequences**.

The principal hotspot positions examined are:

> **R175 · G245 · R248 · R249 · R273 · R282**

The research combines biological sequence resources, similarity-based sequence
identification, protein sequence alignment, residue-level hotspot mapping,
conservation analysis, comparison of canonical and TP53-related sequences,
phylogenetic analysis, and exploratory computational prioritization.

The repository also contains **EleProtect v2.0**, an interactive Streamlit
research application developed around selected computational components of
the project.

---

# 🧭 Research Navigation

| Section | Description |
|---|---|
| [🔬 Research Overview](#-research-overview) | Scientific problem and project scope |
| [🧠 Biological Rationale](#-biological-rationale) | Why TP53 and elephants |
| [🎯 Aim & Objectives](#-aim--objectives) | Research goals |
| [🧬 TP53 Hotspot Framework](#-tp53-hotspot-framework) | Human TP53 positions investigated |
| [🔬 Computational Methodology](#-computational-methodology) | Overall analytical pipeline |
| [🗂️ Data Sources](#️-data-sources) | Biological databases and datasets |
| [📊 Research Components](#-research-components) | Major computational analyses |
| [📓 Research Notebooks](#-research-notebooks) | Interactive computational workflows |
| [🐍 Analysis Scripts](#-analysis-scripts) | Reusable code |
| [📊 Research Outputs](#-research-outputs) | Results and figures |
| [💻 EleProtect v2.0](#-eleprotect-v20) | Interactive research application |
| [📁 Repository Structure](#-repository-structure) | Project organization |
| [♻️ Reproducibility](#️-reproducibility) | Reproducible research framework |
| [⚠️ Interpretation & Limitations](#️-interpretation--limitations) | Scientific boundaries |
| [📚 Research Publication](#-research-publication) | Published research |
| [📄 Research Square](#-research-square) | Broader TP53 research |
| [🌍 ECCB 2026](#-eccb-2026--poster-presentation) | Conference achievement |
| [🎓 Academic Context](#-academic-context) | University and department |
| [🏆 Researcher Achievement](#-researcher-achievement) | Research achievements |
| [🙏 Acknowledgements](#-acknowledgements) | Authors and institutional acknowledgement |
| [👩‍🔬 About the Researcher](#-about-the-researcher) | Researcher profile |
| [🔗 Researcher Profiles](#-researcher-profiles--contact) | LinkedIn, ResearchGate, GitHub and email |
| [📖 Citation](#-citation) | Citation information |
| [📜 License](#-license) | Repository license |
| [🧭 Start Here](#-start-here) | Recommended repository navigation |

---

# 🔬 Research Overview

## The Scientific Question

The central research question is:

> **How conserved are recurrent human TP53 mutation-associated hotspots across
> elephant TP53-related sequences, and what sequence-level patterns emerge
> when canonical elephant TP53 and TP53-related sequences are compared?**

The project connects:

```text
Human Cancer Genomics
        │
        ▼
TP53 Mutation Hotspots
        │
        ▼
Comparative Genomics
        │
        ▼
Elephant TP53-Related Sequences
        │
        ▼
Evolutionary Conservation
        │
        ▼
Cancer Evolution Biology
````

The study is designed as a **comparative and hypothesis-generating
computational analysis**.

The purpose is to identify evolutionary sequence-level patterns that may
support future structural, functional, and experimental investigation.

The research does **not** assume that TP53 sequence conservation alone
explains elephant cancer resistance or provides a complete explanation of
Peto's paradox.

---

# 🧠 Biological Rationale

## Why TP53?

**TP53** encodes the tumour-suppressor protein p53, a major regulator of
cellular responses to genomic stress.

TP53 is involved in:

* DNA-damage response
* Cell-cycle regulation
* Apoptosis
* Cellular senescence
* Genomic stability
* Transcriptional regulation

Recurrent TP53 mutations occur at specific amino-acid positions in human
cancers. These positions provide defined coordinates for residue-level
cross-species comparison.

---

## Why Elephants?

Elephants represent an important comparative system in evolutionary cancer
biology because they combine:

```text
Large body size
       +
Long lifespan
       +
Large number of cells
       ↓
A major evolutionary cancer-biology question
```

This relationship is commonly discussed in the context of **Peto's paradox**.

The project therefore investigates elephant TP53-related sequences as a
comparative system for examining conservation and variation at
human cancer-associated TP53 positions.

The goal is not to reduce elephant cancer biology to a single gene, but to
identify sequence-level patterns that can contribute to broader evolutionary
cancer research.

---

# 🎯 Aim & Objectives

## Aim

To investigate the evolutionary conservation and sequence-level variation of
major human TP53 mutation-associated hotspots across elephant TP53-related
sequences using comparative in-silico bioinformatics approaches.

## Objectives

* Identify recurrent human TP53 mutation-associated hotspots.
* Retrieve human and elephant TP53-related sequence resources.
* Identify TP53-related elephant sequences using similarity-based analysis.
* Curate sequences for comparative analysis.
* Perform protein-level multiple sequence alignment.
* Map human TP53 hotspots onto elephant sequences.
* Evaluate conservation and amino-acid variation at corresponding positions.
* Compare canonical elephant TP53 with TP53-related sequences.
* Examine evolutionary relationships among analysed sequences.
* Develop an exploratory computational prioritization framework.
* Investigate sequence-derived computational and machine-learning approaches.
* Preserve the analysis as a transparent and reproducible research workflow.
* Develop an interactive research interface through **EleProtect v2.0**.

---

# 🧬 TP53 Hotspot Framework

The comparative framework focuses on six major human TP53
mutation-associated positions:

| Hotspot  | Human TP53 Position | Comparative Purpose                    |
| -------- | ------------------: | -------------------------------------- |
| **R175** |                 175 | Mutation-associated structural hotspot |
| **G245** |                 245 | Mutation-associated hotspot            |
| **R248** |                 248 | DNA-binding hotspot                    |
| **R249** |                 249 | Mutation-associated hotspot            |
| **R273** |                 273 | DNA-binding hotspot                    |
| **R282** |                 282 | Structural hotspot                     |

These positions provide defined reference coordinates for residue-level
comparative mapping following sequence alignment.

---

# 🔬 Computational Methodology

The overall computational workflow follows:

```text
Human TP53 Reference
        │
        ▼
Human TP53 Mutation Hotspots
        │
        ▼
Elephant Protein Resources
        │
        ▼
Similarity-Based Sequence Identification
        │
        ▼
Candidate TP53-Related Sequences
        │
        ▼
Sequence Curation
        │
        ▼
Multiple Sequence Alignment
        │
        ▼
Hotspot Position Mapping
        │
        ▼
Conservation / Variation Analysis
        │
        ▼
Canonical vs TP53-Related Comparison
        │
        ▼
Phylogenetic Analysis
        │
        ▼
Computational Prioritization
        │
        ▼
Exploratory Computational / ML Analysis
        │
        ▼
Research Interpretation
```

The detailed methodological description is maintained separately so that
the main README remains focused on the research overview.

👉 **Detailed methodology:**
[`docs/methodology.md`](docs/methodology.md)

👉 **Data provenance:**
[`docs/provenance.md`](docs/provenance.md)

👉 **Reproducibility:**
[`docs/reproducibility.md`](docs/reproducibility.md)

👉 **Interpretation and limitations:**
[`docs/interpretation_and_limitations.md`](docs/interpretation_and_limitations.md)

👉 **Research proposal:**
[`docs/Research_Proposal.pdf`](docs/Research_Proposal.pdf)

---

# 🗂️ Data Sources

The project uses publicly available biological sequence resources.

## Human TP53

| Field             | Information                 |
| ----------------- | --------------------------- |
| Species           | *Homo sapiens*              |
| Gene              | **TP53**                    |
| UniProt accession | **P04637**                  |
| Sequence type     | Canonical protein           |
| Local resource    | `data/raw/human_tp53.fasta` |
| Source            | UniProt                     |

🔗 **UniProt TP53 P04637**
[https://www.uniprot.org/uniprotkb/P04637/entry](https://www.uniprot.org/uniprotkb/P04637/entry)

---

## Asian Elephant Genomic Resource

| Field           | Information                       |
| --------------- | --------------------------------- |
| Species         | *Elephas maximus*                 |
| Assembly        | **GCA_024166365.1**               |
| Assembly        | mEleMax1 primary haplotype        |
| Source          | NCBI Assembly                     |
| Analytical role | Elephant genomic/protein resource |

🔗 **NCBI Assembly**
[https://www.ncbi.nlm.nih.gov/assembly/13211691](https://www.ncbi.nlm.nih.gov/assembly/13211691)

---

## Elephant Proteome

The project uses a large elephant protein FASTA resource for
similarity-based identification of TP53-related sequences.

The complete elephant proteome is approximately **41 MB** and is therefore
not duplicated unnecessarily throughout the repository.

The source, accession, download information and analytical role are preserved
in the project provenance documentation.

👉 [`docs/provenance.md`](docs/provenance.md)

---

# 📊 Research Components

## 1. 🔎 Similarity-Based Sequence Identification

Human TP53 is used as a reference for identifying TP53-related sequences in
the elephant protein resource.

The conceptual workflow is:

```text
Human TP53
    ↓
Similarity Search
    ↓
Elephant Proteome
    ↓
Candidate TP53-Related Proteins
    ↓
Sequence Inspection
    ↓
Curated Comparative Dataset
```

Relevant resources:

👉 [`data/raw/Blast/`](data/raw/Blast/)

---

## 2. 🧬 Multiple Sequence Alignment

Protein-level multiple sequence alignment establishes residue-level
correspondence between human TP53 and selected elephant TP53-related
sequences.

Tools used include:

* **MAFFT**
* **Jalview**

Alignment outputs are maintained under:

👉 [`results/MSA/`](results/MSA/)

👉 [`results/MSA_1and2/`](results/MSA_1and2/)

👉 [`results/MSA_3/`](results/MSA_3/)

---

## 3. 🎯 TP53 Hotspot Mapping

Human TP53 mutation-associated hotspot positions are mapped onto aligned
elephant sequences.

The conceptual process is:

```text
Human TP53 hotspot
        ↓
Reference position
        ↓
Aligned position
        ↓
Corresponding elephant residue
        ↓
Conservation / substitution
        ↓
Comparative interpretation
```

---

## 4. 📊 Conservation Analysis

The project evaluates whether hotspot-associated residues are conserved or
divergent across the analysed sequences.

The interpretation follows:

```text
Observed conservation
        ↓
Evidence consistent with evolutionary constraint
        ↓
Biological hypothesis
```

Observed sequence conservation is not treated as direct proof of biological
function.

---

## 5. 🐘 Canonical TP53 vs TP53-Related Sequences

A component of the project compares canonical elephant TP53 with
TP53-related retrogene-derived sequences.

This enables investigation of sequence-level variation and hotspot-associated
differences among TP53-related sequence categories.

---

## 6. 🌳 Phylogenetic Analysis

Phylogenetic analysis provides evolutionary context for the analysed
sequences.

Relevant outputs:

👉 [`results/phylogeny/`](results/phylogeny/)

👉 [`results/MEGA/`](results/MEGA/)

---

## 7. 📈 Computational Prioritization

The project includes exploratory feature-based computational prioritization.

The framework integrates sequence-derived observations such as similarity,
conservation and mutation relevance to organize candidates for further
investigation.

This prioritization is **research-oriented and exploratory** and should not
be interpreted as a clinical score or validated biological prediction.

👉 [`notebooks/EleProtect_Feature_Prioritization.ipynb`](notebooks/EleProtect_Feature_Prioritization.ipynb)

---

## 8. 🤖 Exploratory Machine Learning

The project contains exploratory computational and machine-learning analyses
based on sequence-derived features.

👉 [`results/ML/`](results/ML/)

👉 [`notebooks/TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb`](notebooks/TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb)

These analyses are exploratory and should be interpreted in the context of
the available dataset, feature representation, validation strategy and
external validation requirements.

---

# 📓 Research Notebooks

The major computational notebooks are:

## `TP53_Comparative_Analysis.ipynb`

Primary comparative TP53 sequence and hotspot analysis.

👉 [`notebooks/TP53_Comparative_Analysis.ipynb`](notebooks/TP53_Comparative_Analysis.ipynb)

---

## `EleProtect_Feature_Prioritization.ipynb`

Feature extraction and exploratory computational prioritization.

👉 [`notebooks/EleProtect_Feature_Prioritization.ipynb`](notebooks/EleProtect_Feature_Prioritization.ipynb)

---

## `TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb`

Exploratory predictive modelling and deep-learning workflow.

👉 [`notebooks/TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb`](notebooks/TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb)

---

# 🐍 Analysis Scripts

Reusable computational implementations are maintained separately from the
interactive notebooks.

👉 [`scripts/`](scripts/)

The principal comparative analysis implementation is:

👉 [`scripts/TP53_Comparative_Analysis.py`](scripts/TP53_Comparative_Analysis.py)

The separation follows:

```text
scripts/
    ↓
Reusable computational implementation

notebooks/
    ↓
Interactive analysis, exploration and interpretation
```

---

# 📊 Research Outputs

The repository preserves computational outputs generated during the
analysis.

## Multiple Sequence Alignments

👉 [`results/MSA/`](results/MSA/)

👉 [`results/MSA_1and2/`](results/MSA_1and2/)

👉 [`results/MSA_3/`](results/MSA_3/)

## Phylogenetic Analysis

👉 [`results/phylogeny/`](results/phylogeny/)

👉 [`results/MEGA/`](results/MEGA/)

## Machine-Learning / Computational Outputs

👉 [`results/ML/`](results/ML/)

## Research Figures

👉 [`figures/`](figures/)

The figures provide visual representations of selected computational
analyses, sequence relationships, alignments, hotspot mapping and
prioritization outputs.

---

# 💻 EleProtect v2.0

## Interactive Research Application

**EleProtect v2.0** is the research-software component associated with this
project.

It provides a Streamlit-based interface for selected sequence-analysis and
exploratory computational-prioritization workflows.

### Features

* 🧬 DNA or protein sequence input
* 🔎 Sequence processing
* 🎯 TP53-oriented analysis
* 📊 Feature extraction
* 🤖 Exploratory ML scoring
* 📁 CSV export
* 🌐 Interactive Streamlit interface

### Application Directory

👉 [`EleProtect_App/`](EleProtect_App/)

### Application Documentation

👉 [`EleProtect_App/README.md`](EleProtect_App/README.md)

### Local Execution

```bash
cd EleProtect_App
pip install -r requirements.txt
streamlit run app.py
```

EleProtect is a **research prototype** and is not intended to function as a
clinical diagnostic, prognostic or therapeutic decision-support system.

---

# 📁 Repository Structure

```text
Elephant-TP53-Hotspot-Mapping/
│
├── data/
│   ├── raw/
│   ├── processed/
│   └── Database/
│
├── notebooks/
│   ├── TP53_Comparative_Analysis.ipynb
│   ├── EleProtect_Feature_Prioritization.ipynb
│   └── TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb
│
├── scripts/
│
├── results/
│   ├── MSA/
│   ├── MSA_1and2/
│   ├── MSA_3/
│   ├── MEGA/
│   ├── ML/
│   └── phylogeny/
│
├── figures/
│
├── docs/
│   ├── README.md
│   ├── methodology.md
│   ├── provenance.md
│   ├── reproducibility.md
│   ├── interpretation_and_limitations.md
│   └── Research_Proposal.pdf
│
├── manuscript/
│
├── EleProtect_App/
│
├── CITATION.cff
├── LICENSE
├── requirements.txt
├── .gitignore
└── README.md
```

---

# 🧭 Where to Find Each Part of the Research

| Interested in                      | Repository location                                                                                                                                    |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 🧬 Raw biological data             | [`data/raw/`](data/raw/)                                                                                                                               |
| 🗂️ Processed datasets             | [`data/processed/`](data/processed/)                                                                                                                   |
| 🔎 BLAST / sequence identification | [`data/raw/Blast/`](data/raw/Blast/)                                                                                                                   |
| 🧬 Comparative TP53 analysis       | [`notebooks/TP53_Comparative_Analysis.ipynb`](notebooks/TP53_Comparative_Analysis.ipynb)                                                               |
| 📈 Feature prioritization          | [`notebooks/EleProtect_Feature_Prioritization.ipynb`](notebooks/EleProtect_Feature_Prioritization.ipynb)                                               |
| 🤖 Predictive modelling            | [`notebooks/TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb`](notebooks/TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb) |
| 🐍 Reusable code                   | [`scripts/`](scripts/)                                                                                                                                 |
| 🧬 Alignment results               | [`results/MSA/`](results/MSA/)                                                                                                                         |
| 🌳 Phylogenetic results            | [`results/phylogeny/`](results/phylogeny/)                                                                                                             |
| 📊 MEGA analysis                   | [`results/MEGA/`](results/MEGA/)                                                                                                                       |
| 🤖 ML outputs                      | [`results/ML/`](results/ML/)                                                                                                                           |
| 🖼️ Figures                        | [`figures/`](figures/)                                                                                                                                 |
| 📖 Methodology                     | [`docs/methodology.md`](docs/methodology.md)                                                                                                           |
| 🧾 Provenance                      | [`docs/provenance.md`](docs/provenance.md)                                                                                                             |
| ♻️ Reproducibility                 | [`docs/reproducibility.md`](docs/reproducibility.md)                                                                                                   |
| ⚠️ Limitations                     | [`docs/interpretation_and_limitations.md`](docs/interpretation_and_limitations.md)                                                                     |
| 📋 Research proposal               | [`docs/Research_Proposal.pdf`](docs/Research_Proposal.pdf)                                                                                             |
| 💻 EleProtect                      | [`EleProtect_App/`](EleProtect_App/)                                                                                                                   |

---

# ♻️ Reproducibility

The repository follows a traceable computational research structure:

```text
Biological Source
       ↓
Data Provenance
       ↓
Raw Data
       ↓
Processing
       ↓
Computational Analysis
       ↓
Results
       ↓
Figures
       ↓
Scientific Interpretation
```

The detailed reproducibility framework is available at:

👉 [`docs/reproducibility.md`](docs/reproducibility.md)

The detailed methodology is available at:

👉 [`docs/methodology.md`](docs/methodology.md)

The data provenance record is available at:

👉 [`docs/provenance.md`](docs/provenance.md)

---

# ⚠️ Interpretation & Limitations

This repository contains computational research and should be interpreted
within the scope of the available sequence-level evidence.

The analysis can provide evidence concerning:

* sequence conservation;
* amino-acid variation;
* sequence similarity;
* TP53-related sequence relationships;
* phylogenetic context; and
* computational prioritization.

The computational analysis does **not independently establish**:

* the complete molecular basis of elephant cancer resistance;
* a complete mechanistic explanation of Peto's paradox;
* functional activity of individual TP53-related retrogenes;
* clinical cancer risk;
* therapeutic response;
* diagnostic utility; or
* experimentally validated biological mechanisms.

The appropriate research progression is:

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

Detailed discussion is available at:

👉 [`docs/interpretation_and_limitations.md`](docs/interpretation_and_limitations.md)

---

# 📚 Research Publication

## *Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants: A Responsible Bioinformatics Innovation Contributing to Cancer Research*

### Authors

**Ritika Rajendra Rawat**
**Sermarani Nadar**
**Gursimran Kaur Uppal**

### Academic Affiliation

**Department of Bioinformatics**
**Guru Nanak Khalsa College**
**University of Mumbai**

### Publication

**2026**

### DOI

[https://doi.org/10.25215/9141002199](https://doi.org/10.25215/9141002199)

### ResearchGate Publication

[https://www.researchgate.net/publication/403581357_COMPARATIVE_IN-SILICO_MAPPING_OF_TP53_MUTATION_HOTSPOTS_IN_ELEPHANTS_A_RESPONSIBLE_BIOINFORMATICS_INNOVATION_CONTRIBUTING_TO_CANCER_RESEARCH](https://www.researchgate.net/publication/403581357_COMPARATIVE_IN-SILICO_MAPPING_OF_TP53_MUTATION_HOTSPOTS_IN_ELEPHANTS_A_RESPONSIBLE_BIOINFORMATICS_INNOVATION_CONTRIBUTING_TO_CANCER_RESEARCH)

This publication represents the focused research direction underlying the
human–elephant TP53 comparative analysis documented in this repository.

---

# 📄 Research Square

## *Evolutionary Conservation and Functional Constraint of TP53 Mutation Hotspots Across Mammalian Species*

The research direction has been extended from the focused human–elephant
comparison toward a broader **mammalian comparative TP53 framework**.

The broader research investigates evolutionary conservation and functional
constraint surrounding recurrent TP53 mutation-associated hotspots across
mammalian species.

### Research Square

[https://doi.org/10.21203/rs.3.rs-9299199/v1](https://doi.org/10.21203/rs.3.rs-9299199/v1)

### ResearchGate

[https://www.researchgate.net/profile/Ritika-Rawat-10](https://www.researchgate.net/profile/Ritika-Rawat-10)

This work represents an extension of the research trajectory from:

```text
Human TP53 ↔ Elephant TP53
```

toward:

```text
TP53 Hotspot Conservation Across Mammalian Species
```

---

# 🌍 ECCB 2026 — Poster Presentation

## European Conference on Computational Biology

The broader TP53 research programme has been accepted for **poster
presentation at ECCB 2026**.

### Poster Title

> **Evolutionary Conservation and Functional Constraint of TP53 Mutation
> Hotspots Across Mammalian Species**

### Conference

**ECCB 2026 — European Conference on Computational Biology**

### Location

**Geneva, Switzerland**

### Dates

**31 August – 4 September 2026**

### Presentation

**Poster Presentation**

🔗 **ECCB 2026**

[https://eccb2026.org/](https://eccb2026.org/)

### Connection to This Research

The ECCB research represents an extension of the research direction documented
in this repository.

The progression is:

```text
Human TP53 Mutation Hotspots
          ↓
Human–Elephant Comparative Analysis
          ↓
TP53 Conservation & Variation
          ↓
Canonical / TP53-Related Sequence Comparison
          ↓
Published Research
          ↓
Broader Mammalian TP53 Analysis
          ↓
Research Square
          ↓
ECCB 2026 Poster Presentation
```

The ECCB presentation therefore represents the broader continuation of the
research into **evolutionary conservation and functional constraint of TP53
mutation-associated hotspots across mammals**.

---

# 🏆 Researcher Achievement

## Ritika Rajendra Rawat

**MSc Bioinformatics**

The project represents a postgraduate computational research direction
connecting comparative genomics, cancer biology and evolutionary
bioinformatics.

### 🏆 ECCB 2026 Poster Presentation

The broader TP53 research:

> **Evolutionary Conservation and Functional Constraint of TP53 Mutation
> Hotspots Across Mammalian Species**

has been accepted for poster presentation at the **European Conference on
Computational Biology (ECCB 2026)** in Geneva, Switzerland.

### 📚 Published Research

Author/co-author of the associated research publication:

> *Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants:
> A Responsible Bioinformatics Innovation Contributing to Cancer Research*

### 📄 Research Square

The research has been extended to a broader mammalian TP53 conservation
framework available through Research Square.

### 💻 Research Software

Developer of **EleProtect v2.0**, an interactive Streamlit-based research
application associated with the TP53 comparative analysis.

### ♻️ Open Computational Research

Maintainer of this GitHub repository containing:

* Research documentation
* Biological data provenance
* Computational notebooks
* Analysis scripts
* Alignment outputs
* Phylogenetic outputs
* Machine-learning outputs
* Research figures
* Reproducibility documentation
* Interactive research software

---

# 🎓 Academic Context

## Department of Bioinformatics

### Guru Nanak Khalsa College

### University of Mumbai

This research was developed within the postgraduate academic and research
context of **Bioinformatics at Guru Nanak Khalsa College, University of
Mumbai**.

The academic environment provided the research and computational biology
context for development of the project.

The associated research publication lists the authors with the Department of
Bioinformatics, Guru Nanak Khalsa College, University of Mumbai.

---

# 🙏 Acknowledgements

The research acknowledges the academic environment and research training
provided by:

### Department of Bioinformatics

**Guru Nanak Khalsa College**
**University of Mumbai**

The project also acknowledges the publicly available biological databases,
sequence resources and computational tools that supported the research,
including:

* **NCBI**
* **UniProt**
* **BLAST**
* **MAFFT**
* **Jalview**
* **MEGA**
* **Python**
* **Google Colab**
* **scikit-learn**
* **Streamlit**

---

## 👥 Research Authors & Contributors

### Ritika Rajendra Rawat

**Primary researcher / author**

Contributed to the computational research framework, comparative analysis,
sequence analysis, research development and associated research outputs.

### Sermarani Nadar

**Co-author**

Co-author of the associated published research on comparative TP53 hotspot
mapping in elephants.

### Gursimran Kaur Uppal

**Co-author**

Co-author of the associated published research on comparative TP53 hotspot
mapping in elephants.

The exact contribution structure should be interpreted according to the
authorship and contribution statements of the associated publication.

---

# 👩‍🔬 About the Researcher

## Ritika Rajendra Rawat

**MSc Bioinformatics**

### Research Interests

* 🧬 Computational Genomics
* 🧬 Comparative Genomics
* 🧬 Cancer Bioinformatics
* 🌳 Evolutionary Bioinformatics
* 🧪 Sequence Analysis
* 🧬 Computational Biology
* 🤖 Machine Learning for Biological Data
* 💻 Research Software
* ♻️ Reproducible Computational Research

The research direction focuses on applying computational approaches to
evolutionary questions in cancer biology, particularly sequence conservation,
tumour-suppressor genes and comparative genomics.

---

# 🔗 Researcher Profiles & Contact

## 💼 LinkedIn

**Ritika Rawat**

[https://in.linkedin.com/in/ritika-rawat-551107219](https://in.linkedin.com/in/ritika-rawat-551107219)

---

## 🔬 ResearchGate

**Ritika Rawat**

[https://www.researchgate.net/profile/Ritika-Rawat-10](https://www.researchgate.net/profile/Ritika-Rawat-10)

---

## 💻 GitHub

**Rita1791**

[https://github.com/Rita1791](https://github.com/Rita1791)

---

## 📧 Research / Academic Email

**[ritikarvl2627@gmail.com](mailto:ritikarvl2627@gmail.com)**

For academic discussions, research collaboration and enquiries regarding this
project.

---

# 📖 Citation

If you use this repository, its methodology, computational framework, figures,
research software or associated research framework, please cite the associated
published work.

### Published Research

```text
Rawat, R. R., Nadar, S., & Uppal, G. K. (2026).

Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants:
A Responsible Bioinformatics Innovation Contributing to Cancer Research.

Ideathon on Vikshit Bharat: Ideas, Innovation and Impact.
Chapter 11, pp. 123–138.

DOI: 10.25215/9141002199
```

Machine-readable citation information is available in:

👉 [`CITATION.cff`](CITATION.cff)

---

# 📜 License

This repository is released under the **MIT License**.

👉 [`LICENSE`](LICENSE)

---

# 🌐 External Research Links

| Resource                    | Link                                                                                                                                                                                                                                                                                                                                                                                               |
| --------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 💻 GitHub Repository        | [https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping](https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping)                                                                                                                                                                                                                                                                             |
| 📚 Published Research       | [https://doi.org/10.25215/9141002199](https://doi.org/10.25215/9141002199)                                                                                                                                                                                                                                                                                                                         |
| 🔬 ResearchGate Publication | [https://www.researchgate.net/publication/403581357_COMPARATIVE_IN-SILICO_MAPPING_OF_TP53_MUTATION_HOTSPOTS_IN_ELEPHANTS_A_RESPONSIBLE_BIOINFORMATICS_INNOVATION_CONTRIBUTING_TO_CANCER_RESEARCH](https://www.researchgate.net/publication/403581357_COMPARATIVE_IN-SILICO_MAPPING_OF_TP53_MUTATION_HOTSPOTS_IN_ELEPHANTS_A_RESPONSIBLE_BIOINFORMATICS_INNOVATION_CONTRIBUTING_TO_CANCER_RESEARCH) |
| 📄 Research Square          | [https://doi.org/10.21203/rs.3.rs-9299199/v1](https://doi.org/10.21203/rs.3.rs-9299199/v1)                                                                                                                                                                                                                                                                                                         |
| 👩‍🔬 ResearchGate Profile  | [https://www.researchgate.net/profile/Ritika-Rawat-10](https://www.researchgate.net/profile/Ritika-Rawat-10)                                                                                                                                                                                                                                                                                       |
| 💼 LinkedIn                 | [https://in.linkedin.com/in/ritika-rawat-551107219](https://in.linkedin.com/in/ritika-rawat-551107219)                                                                                                                                                                                                                                                                                             |
| 🌍 ECCB 2026                | [https://eccb2026.org/](https://eccb2026.org/)                                                                                                                                                                                                                                                                                                                                                     |
| 🧬 UniProt TP53             | [https://www.uniprot.org/uniprotkb/P04637/entry](https://www.uniprot.org/uniprotkb/P04637/entry)                                                                                                                                                                                                                                                                                                   |
| 🐘 NCBI Elephant Assembly   | [https://www.ncbi.nlm.nih.gov/assembly/13211691](https://www.ncbi.nlm.nih.gov/assembly/13211691)                                                                                                                                                                                                                                                                                                   |

---

# 🧭 Start Here

If you are visiting this repository for the first time:

### 1️⃣ Understand the Research

Read this README to understand the scientific question, rationale and
overall research framework.

### 2️⃣ Understand the Methodology

👉 [`docs/methodology.md`](docs/methodology.md)

### 3️⃣ Verify Data Provenance

👉 [`docs/provenance.md`](docs/provenance.md)

### 4️⃣ Inspect the Main Comparative Analysis

👉 [`notebooks/TP53_Comparative_Analysis.ipynb`](notebooks/TP53_Comparative_Analysis.ipynb)

### 5️⃣ Examine Feature Prioritization

👉 [`notebooks/EleProtect_Feature_Prioritization.ipynb`](notebooks/EleProtect_Feature_Prioritization.ipynb)

### 6️⃣ Examine Predictive Modelling

👉 [`notebooks/TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb`](notebooks/TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb)

### 7️⃣ Inspect the Results

👉 [`results/`](results/)

### 8️⃣ Inspect the Figures

👉 [`figures/`](figures/)

### 9️⃣ Review Reproducibility

👉 [`docs/reproducibility.md`](docs/reproducibility.md)

### 🔟 Review Interpretation and Limitations

👉 [`docs/interpretation_and_limitations.md`](docs/interpretation_and_limitations.md)

### 1️⃣1️⃣ Explore EleProtect

👉 [`EleProtect_App/`](EleProtect_App/)

### 1️⃣2️⃣ Explore the Associated Research

📚 **Published Research:**
[https://doi.org/10.25215/9141002199](https://doi.org/10.25215/9141002199)

📄 **Research Square:**
[https://doi.org/10.21203/rs.3.rs-9299199/v1](https://doi.org/10.21203/rs.3.rs-9299199/v1)

🌍 **ECCB 2026:**
[https://eccb2026.org/](https://eccb2026.org/)

---

# 🧬 Research Perspective

This project brings together:

```text
Cancer Biology
      +
Comparative Genomics
      +
Evolutionary Biology
      +
Bioinformatics
      +
Sequence Analysis
      +
Computational Modelling
      +
Research Software
      +
Reproducible Research
```

The central research principle is:

> **Use computational evidence to identify evolutionary patterns, preserve
> the analytical path behind those observations, and translate computational
> findings into hypotheses that can ultimately be investigated through
> structural, functional and experimental research.**

---

<p align="center">

## 🧬 Human TP53 ↔ Elephant TP53 🐘

**Comparative Genomics · Evolutionary Cancer Biology · Computational Biology**

<br><br>

📚 **Published Research**   ·  
📄 **Research Square**   ·  
🌍 **ECCB 2026 Poster Presentation**   ·  
💻 **EleProtect v2.0**

<br><br>

**Ritika Rajendra Rawat · MSc Bioinformatics**

<br><br>

<a href="https://github.com/Rita1791">GitHub</a> · <a href="https://in.linkedin.com/in/ritika-rawat-551107219">LinkedIn</a> · <a href="https://www.researchgate.net/profile/Ritika-Rawat-10">ResearchGate</a>

</p>
```
