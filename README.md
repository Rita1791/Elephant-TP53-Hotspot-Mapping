# 🧬 Elephant TP53 Hotspot Mapping

## Comparative In-Silico Analysis of TP53 Mutation Hotspots Between Humans and Elephants

<p align="center">

### 🧬 Human TP53 &nbsp; ↔ &nbsp; 🐘 Elephant TP53

**Comparative Genomics · Cancer Biology · Evolutionary Bioinformatics · Computational Biology**

</p>

<p align="center">

[![Research](https://img.shields.io/badge/Research-Comparative%20Genomics-1f6feb)](https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping)
[![Python](https://img.shields.io/badge/Python-3.x-3776AB)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/EleProtect-Streamlit-FF4B4B)](https://streamlit.io/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![ECCB 2026](https://img.shields.io/badge/ECCB%202026-Poster%20Presentation-orange)](https://eccb2026.org/)

</p>

---

# 🔬 Research at a Glance

**Elephant TP53 Hotspot Mapping** is an open computational research repository investigating the **evolutionary conservation and sequence-level variation of major human TP53 mutation-associated hotspots across elephant TP53-related sequences**.

The project combines:

- 🧬 Comparative genomics
- 🧪 Cancer biology
- 🌳 Evolutionary bioinformatics
- 🔎 Protein sequence analysis
- 🧬 Multiple sequence alignment
- 🎯 Mutation-hotspot mapping
- 🌳 Phylogenetic analysis
- 📊 Sequence-derived feature analysis
- 🤖 Exploratory machine learning
- 💻 Research software development
- ♻️ Reproducible computational research

The research is motivated by the evolutionary cancer-biology question associated with **Peto's paradox**: why cancer incidence does not simply scale with body size, cell number and lifespan across species.

Rather than assuming that TP53 alone explains elephant cancer resistance, this project uses **comparative sequence analysis to investigate evolutionary patterns that may contribute to future hypotheses about tumour-suppressor evolution**.

---

# 🎯 The Scientific Question

> **How conserved are recurrent human TP53 mutation-associated hotspots across elephant TP53-related sequences, and what sequence-level patterns emerge when canonical elephant TP53 and TP53-related sequences are compared?**

The principal human TP53 hotspot positions investigated are:

| Hotspot | Human TP53 Position | Research Context |
|---|---:|---|
| **R175** | 175 | Structural / mutation-associated hotspot |
| **G245** | 245 | Mutation-associated hotspot |
| **R248** | 248 | DNA-binding hotspot |
| **R249** | 249 | Mutation-associated hotspot |
| **R273** | 273 | DNA-binding hotspot |
| **R282** | 282 | Structural hotspot |

These positions provide defined reference coordinates for cross-species residue-level comparison following sequence alignment.

---

# 🧠 Why This Research?

## Why TP53?

**TP53** encodes the tumour-suppressor protein p53, a central regulator of cellular responses to genomic stress.

TP53 participates in:

- DNA-damage response
- Cell-cycle regulation
- Apoptosis
- Cellular senescence
- Genomic stability
- Transcriptional regulation

Recurrent TP53 mutations occur at specific amino-acid positions in human cancers. These mutation-associated positions provide a defined framework for investigating evolutionary conservation and divergence across species.

---

## 🐘 Why Elephants?

Elephants provide an important comparative system in evolutionary cancer biology because they combine:

```text
Large body size
      +
Long lifespan
      +
Large number of cells
      ↓
A major evolutionary cancer-biology question

This relationship is commonly discussed in the context of Peto's paradox.

The objective of this project is not to claim that elephant cancer resistance is explained by TP53 sequence conservation alone.

Instead, the study asks whether sequence-level patterns surrounding human cancer-associated TP53 hotspots can provide useful evolutionary observations and hypotheses for subsequent structural, functional and experimental investigation.

🧬 Research Framework

The project follows the overall research logic:

Human TP53
    │
    ▼
Human Cancer-Associated Hotspots
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
Conservation & Amino-Acid Variation
    │
    ▼
Canonical vs TP53-Related Comparison
    │
    ▼
Phylogenetic Analysis
    │
    ▼
Sequence-Derived Feature Analysis
    │
    ▼
Exploratory Computational Prioritization
    │
    ▼
Biological Interpretation
🎯 Aim & Objectives
Aim

To investigate the evolutionary conservation and sequence-level variation of major human TP53 mutation-associated hotspots across elephant TP53-related sequences using comparative in-silico bioinformatics approaches.

Objectives
Identify recurrent human TP53 mutation-associated hotspots.
Retrieve human and elephant TP53-related sequence resources.
Identify TP53-related elephant sequences using similarity-based analysis.
Curate sequences for comparative analysis.
Perform protein-level multiple sequence alignment.
Map human TP53 hotspots onto elephant sequences.
Evaluate conservation and amino-acid variation.
Compare canonical elephant TP53 with TP53-related sequences.
Examine evolutionary relationships among analysed sequences.
Develop exploratory sequence-derived computational prioritization.
Investigate exploratory machine-learning approaches.
Preserve the computational workflow as a transparent and reproducible research resource.
Develop EleProtect v2.0, an interactive research application associated with the analysis.
📊 What Was Analysed?

The repository contains computational work covering:

1. 🔎 Sequence Identification

Human TP53 is used as a reference to identify TP53-related sequences within elephant protein resources.

2. 🧬 Multiple Sequence Alignment

Protein sequences are aligned to establish residue-level correspondence between human TP53 and elephant TP53-related sequences.

Primary tools include:

MAFFT
Jalview
3. 🎯 Hotspot Mapping

Human TP53 mutation-associated positions are mapped onto aligned elephant sequences.

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
4. 📈 Conservation Analysis

The analysis examines conservation and amino-acid variation at positions corresponding to human TP53 hotspots.

5. 🐘 TP53-Related Sequence Comparison

Canonical elephant TP53 and TP53-related / retrogene-derived sequences are compared to investigate sequence-level variation.

6. 🌳 Phylogenetic Analysis

Phylogenetic analysis provides evolutionary context for the analysed sequences.

7. 🤖 Exploratory Computational / ML Analysis

Sequence-derived features are used in exploratory computational prioritization and machine-learning workflows.

These analyses are exploratory and are not presented as clinically validated prediction models.

📁 Repository Navigation
Research Component	Location
🧬 Raw biological data	data/raw/
🗂️ Processed data	data/processed/
🔎 Sequence identification	data/raw/Blast/
📓 Main comparative analysis	notebooks/TP53_Comparative_Analysis.ipynb
📈 Feature prioritization	notebooks/EleProtect_Feature_Prioritization.ipynb
🤖 Predictive modelling	notebooks/TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb
🐍 Reusable analysis scripts	scripts/
🧬 MSA results	results/MSA/
🧬 Additional MSA results	results/MSA_1and2/
🧬 Three-sequence MSA	results/MSA_3/
🌳 Phylogeny	results/phylogeny/
📊 MEGA outputs	results/MEGA/
🤖 ML outputs	results/ML/
🖼️ Research figures	figures/
📖 Methodology	docs/methodology.md
🧾 Data provenance	docs/provenance.md
♻️ Reproducibility	docs/reproducibility.md
⚠️ Interpretation & limitations	docs/interpretation_and_limitations.md
📋 Research proposal	docs/Research_Proposal.pdf
💻 EleProtect application	EleProtect_App/
📓 Main Research Notebooks
TP53_Comparative_Analysis.ipynb

Primary computational workflow for the comparative TP53 analysis.

👉 Open notebook

EleProtect_Feature_Prioritization.ipynb

Feature extraction and exploratory computational prioritization.

👉 Open notebook

TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb

Exploratory computational and predictive-modelling workflow.

👉 Open notebook

🐍 Reusable Computational Code

Reusable implementations are maintained under:

👉 scripts/

The main comparative implementation is:

👉 scripts/TP53_Comparative_Analysis.py

The repository separates:

scripts/
    ↓
Reusable computational implementation

notebooks/
    ↓
Interactive analysis, exploration and interpretation

This separation helps distinguish exploratory research from reusable computational components.

📊 Research Outputs

The repository preserves computational outputs generated during the study.

🧬 Multiple Sequence Alignment

results/MSA/

results/MSA_1and2/

results/MSA_3/

🌳 Phylogenetic Analysis

results/phylogeny/

results/MEGA/

🤖 Computational / ML Outputs

results/ML/

🖼️ Figures

figures/

📋 Results Documentation

results/README.md

♻️ Reproducibility

The repository is structured to preserve a traceable relationship between:

Biological Source
      ↓
Data Provenance
      ↓
Raw Data
      ↓
Processing
      ↓
Analysis
      ↓
Results
      ↓
Figures
      ↓
Interpretation

Detailed documentation:

docs/methodology.md
docs/provenance.md
docs/reproducibility.md
docs/interpretation_and_limitations.md
⚠️ Scientific Scope & Limitations

This repository contains computational research and should be interpreted within the scope of the available sequence-level evidence.

The analysis can provide evidence concerning:

sequence similarity
sequence conservation
amino-acid variation
TP53-related sequence relationships
phylogenetic context
computational prioritization

The analysis does not independently establish:

the complete molecular basis of elephant cancer resistance;
a complete mechanistic explanation of Peto's paradox;
functional activity of individual TP53-related retrogenes;
clinical cancer risk;
therapeutic response;
diagnostic utility; or
experimentally validated biological mechanisms.

The intended scientific progression is:

Computational Observation
        ↓
Evolutionary Interpretation
        ↓
Hypothesis Generation
        ↓
Structural / Functional Investigation
        ↓
Experimental Validation

For the detailed discussion:

👉 docs/interpretation_and_limitations.md

💻 EleProtect v2.0
Interactive Research Application

EleProtect v2.0 is the software component associated with this research project.

It provides a Streamlit-based interface for selected sequence-analysis and exploratory computational-prioritization workflows.

Features
🧬 DNA / protein sequence input
🔎 Sequence processing
🎯 TP53-oriented analysis
📊 Feature extraction
🤖 Exploratory ML scoring
📁 CSV export
🌐 Interactive Streamlit interface
Application

👉 EleProtect_App/

👉 EleProtect_App/README.md

Local execution
cd EleProtect_App
pip install -r requirements.txt
streamlit run app.py

Important: EleProtect is a research prototype and is not intended to function as a clinical diagnostic, prognostic, or therapeutic decision-support system.

📚 Associated Research Publication
Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants: A Responsible Bioinformatics Innovation Contributing to Cancer Research
Authors

Ritika Rajendra Rawat
Sermarani Nadar
Gursimran Kaur Uppal

Academic Affiliation

Department of Bioinformatics
Guru Nanak Khalsa College of Arts, Commerce and Science
University of Mumbai

Publication

2026

DOI

10.25215/9141002199

ResearchGate

View the research profile

📄 Research Square

The broader TP53 research direction is available through Research Square:

Research Square DOI:
10.21203/rs.3.rs-9299199/v1

The broader research direction extends the human–elephant comparison toward investigation of TP53 hotspot conservation and functional constraint across mammals.

🌍 ECCB 2026 — Poster Presentation
European Conference on Computational Biology

The broader TP53 research programme has been accepted for poster presentation at ECCB 2026.

Poster

Evolutionary Conservation and Functional Constraint of TP53 Mutation Hotspots Across Mammalian Species

Conference

ECCB 2026 — European Conference on Computational Biology

Location

📍 Geneva, Switzerland

Dates

📅 31 August – 4 September 2026

Presentation Format

🎓 Poster Presentation

🔗 ECCB 2026 Conference

🔗 How ECCB Connects to This Repository

The ECCB work represents a broader continuation of the research direction developed through this project.

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

This repository therefore represents one component of a continuing research trajectory in comparative genomics, evolutionary cancer biology and computational analysis of TP53.

🎓 Academic Context

This research was developed within the postgraduate academic and research context of:

Department of Bioinformatics

Guru Nanak Khalsa College of Arts, Commerce and Science

University of Mumbai

The project represents a computational research direction developed through postgraduate training in Bioinformatics, connecting:

bioinformatics → comparative genomics → cancer biology → evolutionary analysis → research software

🏆 Researcher Profile
👩‍🔬 Ritika Rajendra Rawat

MSc Bioinformatics

Computational biology researcher working at the intersection of:

🧬 Computational Genomics
🧬 Comparative Genomics
🧪 Cancer Bioinformatics
🌳 Evolutionary Bioinformatics
🧬 Sequence Analysis
🤖 Machine Learning for Biological Data
💻 Research Software
♻️ Reproducible Computational Research
Research Direction

My research interests focus on applying computational approaches to biological and evolutionary questions, particularly:

How can sequence-level computational evidence reveal evolutionary constraints relevant to cancer biology?

This project represents an effort to build an open, inspectable and reproducible computational research framework rather than treating computational results as isolated findings.

🏅 Research Highlights
Achievement	Details
🎓 Academic	MSc Bioinformatics
📚 Publication	Co-authored published TP53 comparative research
🌍 ECCB 2026	Accepted poster presentation
🧬 Research	Comparative TP53 hotspot analysis
💻 Software	Developer of EleProtect v2.0
♻️ Open Research	Reproducible GitHub research repository
📄 Research Square	Broader mammalian TP53 research
🔗 Researcher & Contact
💼 LinkedIn

Ritika Rawat

👉 linkedin.com/in/ritika-rawat-551107219

🔬 ResearchGate

Ritika Rawat

👉 ResearchGate Profile

💻 GitHub

Rita1791

👉 github.com/Rita1791

📧 Academic / Research Email

ritika.rawat27@outlook.com

For:

academic collaboration
research discussions
computational biology projects
PhD-related academic communication
research enquiries
🙏 Acknowledgements

The project acknowledges the academic and research environment provided by:

Department of Bioinformatics

Guru Nanak Khalsa College of Arts, Commerce and Science
University of Mumbai

The project also acknowledges the publicly available biological resources and computational tools that supported the analysis, including:

NCBI
UniProt
BLAST
MAFFT
Jalview
MEGA
Python
Google Colab
scikit-learn
Streamlit
👥 Authors & Contributors
Ritika Rajendra Rawat

Primary researcher / author

Contributed to the computational research framework, comparative analysis, sequence analysis, research development and associated research outputs.

Sermarani Nadar

Co-author

Co-author of the associated published research on comparative TP53 hotspot mapping in elephants.

Gursimran Kaur Uppal

Co-author

Co-author of the associated published research on comparative TP53 hotspot mapping in elephants.

Authorship and contribution should ultimately be interpreted according to the formal contribution statement associated with the publication.

🗂️ Data Resources
Human TP53
Species: Homo sapiens
Gene: TP53
UniProt accession: P04637
Sequence type: Canonical protein

👉 UniProt P04637 — TP53

Elephant Genomic Resource
Species: Elephas maximus
Assembly: GCA_024166365.1
Assembly name: mEleMax1 primary haplotype
Source: NCBI Assembly

👉 NCBI Assembly

Detailed provenance:

👉 docs/provenance.md

🧭 Recommended Path for a Researcher / PI

If you are reviewing this repository for the first time:

01 — Understand the research

Read this README.

02 — Inspect the methodology

👉 docs/methodology.md

03 — Verify data provenance

👉 docs/provenance.md

04 — Inspect reproducibility

👉 docs/reproducibility.md

05 — Open the main analysis

👉 notebooks/TP53_Comparative_Analysis.ipynb

06 — Inspect computational outputs

👉 results/

07 — Examine figures

👉 figures/

08 — Examine exploratory prioritization

👉 notebooks/EleProtect_Feature_Prioritization.ipynb

09 — Explore the research application

👉 EleProtect_App/

10 — Read the limitations

👉 docs/interpretation_and_limitations.md

11 — Explore associated research

📚 Published Research

📄 Research Square

🌍 ECCB 2026

📖 Citation

If you use this repository, its methodology, computational framework, figures, research software, or associated research framework, please cite the associated publication.

Rawat, R. R., Nadar, S., & Uppal, G. K. (2026).

Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants:
A Responsible Bioinformatics Innovation Contributing to Cancer Research.

Ideathon on Vikshit Bharat: Ideas, Innovation and Impact.
Chapter 11, pp. 123–138.

DOI: 10.25215/9141002199

Machine-readable citation information:

👉 CITATION.cff

📜 License

This repository is released under the MIT License.

👉 LICENSE

🔗 External Research Links
Resource	Link
💻 GitHub Repository	Elephant-TP53-Hotspot-Mapping
📚 Published Research	DOI: 10.25215/9141002199
📄 Research Square	DOI: 10.21203/rs.3.rs-9299199/v1
🔬 ResearchGate Profile	Ritika Rawat
💼 LinkedIn	Ritika Rawat
📧 Email	ritika.rawat27@outlook.com
🌍 ECCB 2026	European Conference on Computational Biology
🧬 UniProt TP53	P04637
🐘 NCBI Elephant Assembly	GCA_024166365.1
🌐 Research Perspective

This project brings together:

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

The central principle of this repository is:

Use computational evidence to identify evolutionary patterns, preserve the analytical path behind those observations, and translate computational findings into hypotheses that can ultimately be investigated through structural, functional and experimental research.

<p align="center">
🧬 Human TP53 ↔ Elephant TP53 🐘

Comparative Genomics · Evolutionary Cancer Biology · Computational Biology

<br>

📚 Published Research   ·  
📄 Research Square   ·  
🌍 ECCB 2026 Poster Presentation   ·  
💻 EleProtect v2.0

<br><br>

Ritika Rajendra Rawat · MSc Bioinformatics

<br><br>

<a href="https://github.com/Rita1791">GitHub</a> ·
<a href="https://in.linkedin.com/in/ritika-rawat-551107219">LinkedIn</a> ·
<a href="https://www.researchgate.net/profile/Ritika-Rawat-10">ResearchGate</a> ·
<a href="mailto:ritika.rawat27@outlook.com">Email</a>

</p> ```
