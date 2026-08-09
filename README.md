# 🧬 Elephant TP53 Hotspot Mapping

## 🧑 Human TP53 ↔ 🐘 Elephant TP53

### Comparative In-Silico Analysis of TP53 Mutation Hotspots Between Humans and Elephants

**Comparative Genomics · Cancer Biology · Evolutionary Bioinformatics · Computational Biology · Reproducible Research**

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

The project integrates:

- 🧬 Comparative genomics
- 🧪 Cancer biology
- 🌳 Evolutionary bioinformatics
- 🧑 Human TP53 reference analysis
- 🐘 Elephant TP53 sequence analysis
- 🔎 Sequence identification
- 🧬 Multiple sequence alignment
- 🎯 Mutation-hotspot mapping
- 📊 Conservation analysis
- 🌳 Phylogenetic analysis
- 🤖 Exploratory machine-learning approaches
- 💻 Research software development
- ♻️ Reproducible computational research

The project is motivated by the evolutionary cancer-biology question associated with **Peto's paradox**: why cancer incidence does not simply scale with body size, cell number, and lifespan across species.

The study does **not** assume that TP53 sequence conservation alone explains elephant cancer resistance. Instead, it investigates sequence-level evolutionary patterns that can contribute to the generation of hypotheses for subsequent structural, functional, and experimental research.

---

# 🎯 Scientific Question

> **How conserved are recurrent human TP53 mutation-associated hotspots across elephant TP53-related sequences, and what sequence-level patterns emerge when canonical elephant TP53 and TP53-related sequences are compared?**

The principal human TP53 hotspot positions investigated in this project include:

| Hotspot | Human TP53 Position | Research Context |
|---|---:|---|
| **R175** | 175 | Structural / mutation-associated hotspot |
| **G245** | 245 | Mutation-associated hotspot |
| **R248** | 248 | DNA-binding hotspot |
| **R249** | 249 | Mutation-associated hotspot |
| **R273** | 273 | DNA-binding hotspot |
| **R282** | 282 | Structural hotspot |

These positions provide defined human reference coordinates for cross-species residue-level comparison following sequence alignment.

---

# 🧠 Why This Research?

## 🧑 Why Human TP53?

**TP53** encodes the tumour-suppressor protein p53, a central regulator of cellular responses to genomic stress.

TP53 is involved in:

- DNA-damage response
- Cell-cycle regulation
- Apoptosis
- Cellular senescence
- Genomic stability
- Transcriptional regulation

Recurrent TP53 mutations occur at functionally important positions in human cancers. These mutation-associated hotspots provide a defined framework for investigating whether corresponding residues are conserved or divergent across species.

### 🧬 Human TP53 Reference

| Attribute | Details |
|---|---|
| Species | *Homo sapiens* |
| Gene | **TP53** |
| UniProt Accession | **P04637** |
| Sequence Type | Canonical protein |

🔗 **[UniProt P04637 — Human TP53](https://www.uniprot.org/uniprotkb/P04637/entry)**

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

The relationship between body size, lifespan, cell number, and cancer incidence is commonly discussed in the context of Peto's paradox.

This project therefore uses elephants as a comparative evolutionary system to investigate TP53 sequence conservation and variation.

The research does not claim that elephant cancer resistance is explained by TP53 sequence conservation alone.

Instead, the computational objective is to identify sequence-level patterns that may inform future biological investigation.

🧬 Research Framework

The overall analytical logic is:

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
🔬 Computational Methodology

The project follows a multi-stage comparative bioinformatics workflow.

1. Sequence Identification

Human TP53 is used as a reference for identifying TP53-related sequences within elephant protein resources.

Similarity-based analysis is used to identify candidate sequences for subsequent evaluation.

2. Sequence Curation

Candidate sequences are inspected and organized into appropriate datasets for comparative analysis.

The repository distinguishes between:

Raw biological resources
Processed sequences
Alignment-ready datasets
Analytical outputs
3. Multiple Sequence Alignment

Protein sequences are aligned to establish residue-level correspondence between human TP53 and elephant TP53-related sequences.

Primary computational resources include:

MAFFT
Jalview
MEGA
4. Hotspot Mapping

Human TP53 mutation-associated positions are mapped onto aligned sequences.

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
5. Conservation Analysis

The analysis evaluates sequence conservation and amino-acid variation at positions corresponding to selected human TP53 hotspots.

6. Canonical vs TP53-Related Sequence Comparison

Canonical elephant TP53 and TP53-related / retrogene-derived sequences are compared to investigate sequence-level variation.

7. Phylogenetic Analysis

Phylogenetic analysis provides evolutionary context for the analysed sequences.

8. Sequence-Derived Feature Analysis

Sequence-derived characteristics are calculated for analysed sequences, including composition-related features and similarity measures.

9. Exploratory Computational Prioritization

Sequence-derived features are used in exploratory computational and machine-learning workflows.

These analyses are exploratory and are not presented as clinically validated prediction models.

🗂️ Data Sources

The project uses publicly accessible biological sequence resources.

Resource	Purpose
NCBI	Biological sequence and genomic resources
UniProt	Protein sequence and annotation resources
NCBI BLAST	Similarity-based sequence identification
Elephant genomic resources	Comparative elephant sequence analysis
Human TP53 resources	Reference sequence and hotspot analysis
🐘 Elephant Genomic Reference
Attribute	Details
Species	Elephas maximus
Assembly	GCA_024166365.1
Assembly Name	mEleMax1 primary haplotype
Source	NCBI Assembly

🔗 NCBI Assembly — GCA_024166365.1

For complete sequence provenance, accession information, and source documentation:

👉 Read docs/provenance.md

📁 Repository Structure
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
🧭 Researcher / PI Navigation

If you are a researcher, professor, PI, collaborator, or reviewer visiting this repository for the first time, the recommended path is:

01 — 🔬 Understand the Research

Start with this README to understand:

Scientific question
Biological rationale
Research objectives
Computational framework
Research outputs

👉 Start with this README

02 — 🧪 Inspect the Methodology

Review the detailed computational methodology:

👉 docs/methodology.md

03 — 🧾 Verify Data Provenance

Inspect biological sources, accession information, and data provenance:

👉 docs/provenance.md

04 — ♻️ Inspect Reproducibility

Review the computational reproduction framework:

👉 docs/reproducibility.md

05 — 📓 Open the Main Analysis

Explore the principal comparative analysis:

👉 TP53_Comparative_Analysis.ipynb

06 — 📊 Inspect Results

Explore computational outputs:

👉 results/

07 — 🖼️ Examine Figures

Explore research visualizations:

👉 figures/

08 — 🤖 Explore Feature Prioritization

Review the exploratory sequence-feature analysis:

👉 EleProtect_Feature_Prioritization.ipynb

09 — 💻 Explore EleProtect

Explore the associated research software:

👉 EleProtect_App/

10 — ⚠️ Review Scientific Limitations

Review the interpretation boundaries of the computational findings:

👉 docs/interpretation_and_limitations.md

11 — 📚 Explore Associated Research

See the publication, Research Square work, and ECCB research direction below.

📊 Research Outputs

The repository preserves computational outputs generated during the study.

🧬 Multiple Sequence Alignment

👉 results/MSA/

👉 results/MSA_1and2/

👉 results/MSA_3/

🌳 Phylogenetic Analysis

👉 results/phylogeny/

👉 results/MEGA/

🤖 Computational / ML Outputs

👉 results/ML/

🖼️ Research Figures

👉 figures/

📋 Results Documentation

👉 results/README.md

♻️ Reproducibility

A central objective of this repository is to maintain a traceable relationship between the biological source data, computational analysis, and resulting observations.

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

Detailed documentation:

🧪 Methodology
🧾 Data Provenance
♻️ Reproducibility
⚠️ Interpretation & Limitations
🐍 Reusable Computational Code

Reusable computational implementations are maintained under:

👉 scripts/

The principal comparative implementation is:

👉 scripts/TP53_Comparative_Analysis.py

The repository separates reusable computational implementation from interactive analysis:

scripts/
    ↓
Reusable computational implementation

notebooks/
    ↓
Interactive analysis, exploration and interpretation

This separation is intended to make the research workflow easier to inspect, reuse, and extend.

📓 Main Research Notebooks
1. 🧬 TP53 Comparative Analysis

Primary computational workflow for the comparative TP53 analysis.

👉 Open TP53_Comparative_Analysis.ipynb

2. 🤖 EleProtect Feature Prioritization

Feature extraction and exploratory computational prioritization.

👉 Open EleProtect_Feature_Prioritization.ipynb

3. 🧠 Predictive Modelling

Exploratory computational and predictive-modelling workflow.

👉 Open TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb

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
Explore the Application

👉 Open EleProtect_App

👉 Read the EleProtect documentation

Local Execution
cd EleProtect_App
pip install -r requirements.txt
streamlit run app.py

⚠️ Research-use statement: EleProtect is a research prototype and is not intended to function as a clinical diagnostic, prognostic, or therapeutic decision-support system.

⚠️ Scientific Scope & Limitations

This repository contains computational research and should be interpreted within the scope of the available sequence-level evidence.

The analysis can provide evidence concerning:

Sequence similarity
Sequence conservation
Amino-acid variation
TP53-related sequence relationships
Phylogenetic context
Sequence-derived computational prioritization

The analysis does not independently establish:

The complete molecular basis of elephant cancer resistance
A complete mechanistic explanation of Peto's paradox
Functional activity of individual TP53-related retrogenes
Clinical cancer risk
Therapeutic response
Diagnostic utility
Experimentally validated biological mechanisms

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

For the detailed scientific discussion:

👉 Read Interpretation & Limitations

📚 Associated Research Publication
Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants

Rawat, R. R., Nadar, S., & Uppal, G. K. (2026).

Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants: A Responsible Bioinformatics Innovation Contributing to Cancer Research.

Ideathon on Vikshit Bharat: Ideas, Innovation and Impact.

Chapter 11, pp. 123–138.

🔗 Publication

DOI: 10.25215/9141002199

🧾 Citation Metadata

👉 CITATION.cff

If you use this repository, its methodology, computational framework, figures, research software, or associated research framework, please cite the associated publication.

📄 Research Square

The broader TP53 research direction extends this comparative work toward investigation of TP53 hotspot conservation and functional constraint across mammalian species.

🔗 Research Square

DOI: 10.21203/rs.3.rs-9299199/v1

This research direction provides a broader mammalian context for the evolutionary questions explored in the present repository.

🌍 ECCB 2026 — Poster Presentation
European Conference on Computational Biology

The broader TP53 research programme has been accepted for poster presentation at ECCB 2026.

Poster Title

Evolutionary Conservation and Functional Constraint of TP53 Mutation Hotspots Across Mammalian Species

Conference

ECCB 2026 — European Conference on Computational Biology

📍 Geneva, Switzerland

📅 31 August – 4 September 2026

🎓 Poster Presentation

🔗 ECCB 2026 Conference

🔗 Relationship Between This Repository and ECCB

The ECCB work represents a broader continuation of the TP53 research direction developed through this project.

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

The present repository therefore represents one component of a continuing research trajectory in:

Comparative Genomics · Evolutionary Cancer Biology · Computational Biology

🎓 Academic Context

This research was developed within the postgraduate academic and research context of:

Department of Bioinformatics

Guru Nanak Khalsa College of Arts, Commerce and Science

University of Mumbai

The project represents a computational research direction connecting:

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
🙏 Acknowledgements

The project acknowledges the academic environment associated with:

Department of Bioinformatics
Guru Nanak Khalsa College of Arts, Commerce and Science
University of Mumbai

The project also acknowledges the publicly available biological resources and computational tools that supported the analysis:

NCBI
UniProt
NCBI BLAST
MAFFT
Jalview
MEGA
Python
Google Colab
scikit-learn
Streamlit
👥 Authors & Contributors
👩‍🔬 Ritika Rajendra Rawat

Primary Researcher / Author

Contributed to the computational research framework, comparative analysis, sequence analysis, research development, and associated research outputs.

🧑‍🔬 Sermarani Nadar

Co-author

Co-author of the associated published research on comparative TP53 hotspot mapping in elephants.

🧑‍🔬 Gursimran Kaur Uppal

Co-author

Co-author of the associated published research on comparative TP53 hotspot mapping in elephants.

Authorship and contribution should ultimately be interpreted according to the formal contribution statement associated with the publication.

🏆 Researcher Profile
👩‍🔬 Ritika Rajendra Rawat

MSc Bioinformatics

Computational biology researcher working at the intersection of:

🧬 Computational Genomics
🧬 Comparative Genomics
🧪 Cancer Bioinformatics
🌳 Evolutionary Bioinformatics
🧬 Protein & Sequence Analysis
🤖 Machine Learning for Biological Data
💻 Research Software
♻️ Reproducible Computational Research
🔬 Research Direction

My research interests focus on applying computational approaches to biological and evolutionary questions, particularly:

How can sequence-level computational evidence reveal evolutionary constraints relevant to cancer biology?

This project represents an effort to build an open, inspectable, and reproducible computational research framework rather than treating computational results as isolated findings.

🏅 Research Highlights
Achievement	Details
🎓 Academic	MSc Bioinformatics
📚 Publication	Co-authored published TP53 comparative research
🌍 ECCB 2026	Accepted poster presentation
🧬 Research	Comparative TP53 hotspot analysis
💻 Software	Developer of EleProtect v2.0
♻️ Open Research	Reproducible GitHub research repository
📄 Research Square	Broader mammalian TP53 research
🔗 Researcher & Professional Links
💻 GitHub

Rita1791

👉 GitHub Profile

👉 This Research Repository

💼 LinkedIn

Ritika Rajendra Rawat

👉 Connect with me on LinkedIn

For:

Academic networking
Research collaboration
Computational biology discussions
PhD-related academic communication
Scientific networking
🔬 ResearchGate

Ritika Rawat

👉 View my ResearchGate profile

ResearchGate provides access to the researcher's broader research and publication profile.

📧 Academic / Research Email

ritika.rawat27@outlook.com

👉 Contact me by email

Suitable for:

Academic collaboration
Research discussions
Computational biology projects
PhD-related academic communication
Research enquiries
📖 Citation

If you use this repository, its methodology, computational framework, figures, research software, or associated research framework, please cite:

Rawat, R. R., Nadar, S., & Uppal, G. K. (2026).
Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants: A Responsible Bioinformatics Innovation Contributing to Cancer Research.
Ideathon on Vikshit Bharat: Ideas, Innovation and Impact.
Chapter 11, pp. 123–138.
DOI: 10.25215/9141002199

Machine-Readable Citation

👉 CITATION.cff

📜 License

This repository is released under the MIT License.

👉 View LICENSE

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

Use computational evidence to identify evolutionary patterns, preserve the analytical path behind those observations, and translate computational findings into hypotheses that can ultimately be investigated through structural, functional, and experimental research.

🧭 Research Philosophy

The project is designed around a simple principle:

Transparent Data
      +
Traceable Analysis
      +
Reproducible Computation
      +
Responsible Interpretation
      =
Open Computational Research

The repository therefore preserves not only selected results, but also the computational pathway through which those results were generated.

🚀 Future Research Direction

Potential future extensions include:

Structural analysis of conserved and divergent TP53 regions
Comparative analysis across additional mammalian species
Functional investigation of TP53-related sequences
Integration of structural and evolutionary features
Expanded comparative cancer-genomics datasets
More rigorous benchmarking of computational prioritization methods
Experimental validation of biologically relevant hypotheses

The intended progression is:

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
🔗 Quick Access
Research Resource	Open
🧬 Repository	Elephant-TP53-Hotspot-Mapping
🧪 Methodology	docs/methodology.md
🧾 Provenance	docs/provenance.md
♻️ Reproducibility	docs/reproducibility.md
⚠️ Limitations	docs/interpretation_and_limitations.md
📊 Results	results/
🖼️ Figures	figures/
📓 Main Notebook	TP53_Comparative_Analysis.ipynb
🤖 Feature Prioritization	EleProtect_Feature_Prioritization.ipynb
💻 EleProtect	EleProtect_App/
📚 Publication	DOI 10.25215/9141002199
📄 Research Square	DOI 10.21203/rs.3.rs-9299199/v1
🌍 ECCB 2026	eccb2026.org
🔬 ResearchGate	Ritika Rawat
💼 LinkedIn	Ritika Rawat
📧 Email	ritika.rawat27@outlook.com
<p align="center">
🧑 Human TP53 ↔ 🐘 Elephant TP53

Comparative Genomics · Evolutionary Cancer Biology · Computational Biology

<br>

📚 Published Research
  ·  
📄 Research Square
  ·  
🌍 ECCB 2026
  ·  
💻 EleProtect v2.0

<br><br>

👩‍🔬 Ritika Rajendra Rawat

MSc Bioinformatics

<br>

<a href="https://github.com/Rita1791">GitHub</a>
  ·  
<a href="https://in.linkedin.com/in/ritika-rawat-551107219">LinkedIn</a>
  ·  
<a href="https://www.researchgate.net/profile/Ritika-Rawat-10">ResearchGate</a>
  ·  
<a href="mailto:ritika.rawat27@outlook.com">Email</a>

</p> ```
