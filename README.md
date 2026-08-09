# 🧬 Elephant TP53 Hotspot Mapping

## Comparative In-Silico Analysis of TP53 Mutation Hotspots Between Humans and Elephants

<p align="center">

<b>Comparative Genomics • Cancer Biology • Evolutionary Bioinformatics • Reproducible Research Software</b>

</p>

<p align="center">

🧬 Human TP53 &nbsp; ↔ &nbsp; 🐘 Elephant TP53  
<br>
🔎 Hotspot Mapping &nbsp; • &nbsp; 🧪 Sequence Analysis &nbsp; • &nbsp; 🌳 Evolutionary Context  
<br>
🤖 Exploratory Computational Modelling &nbsp; • &nbsp; 💻 EleProtect

</p>

---

## 📑 Research Navigation

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
|---|---|
| 🧬 **Research Area** | Comparative Genomics • Cancer Biology • Evolutionary Bioinformatics |
| 🧬 **Primary Gene** | **TP53 / p53** |
| 🐘 **Comparative System** | Human • Elephant |
| 🎯 **Focal Hotspots** | R175 • G245 • R248 • R249 • R273 • R282 |
| 🔬 **Core Analysis** | Protein sequence comparison • Multiple sequence alignment • Hotspot mapping • Conservation analysis |
| 🧪 **Computational Framework** | Sequence analysis • Similarity analysis • Alignment • Phylogenetics • Feature analysis |
| 🤖 **Exploratory Component** | Sequence-derived feature analysis and machine-learning ranking |
| 📊 **Research Status** | ✅ Completed research project |
| 📚 **Publication Status** | ✅ Published |
| 🐘 **Companion Research Software** | **EleProtect v2.0** |
| 🎓 **Research Context** | Postgraduate Bioinformatics Research |

> **Scientific positioning:** This repository is a computational research
> framework for comparative TP53 sequence analysis and hypothesis generation.
> It is not a clinically validated cancer prediction system.

---

# 📄 Abstract

This research presents a comparative in-silico analysis of **TP53 mutation
hotspots** between humans and elephants, with emphasis on sequence conservation,
variation, and comparative analysis of canonical and TP53-related elephant
sequences.

The study was motivated by the evolutionary question associated with
**Peto's paradox**, which concerns the lack of a simple proportional
relationship between body size, lifespan, cell number, and cancer incidence
across species.

Recurrent human TP53 mutation-associated positions were examined within a
comparative sequence framework to investigate whether corresponding residues
show conservation or variation in elephant TP53-related sequences.

The computational workflow integrates publicly available biological resources,
sequence preprocessing, similarity-based sequence identification, protein-level
sequence alignment, mutation-associated hotspot mapping, conservation analysis,
phylogenetic analysis, sequence-derived feature analysis, and exploratory
computational modelling.

The resulting framework provides a reproducible computational basis for
examining evolutionary patterns surrounding TP53 and for generating hypotheses
that can be investigated through future structural, functional, comparative,
and experimental studies.

---

# 🧠 Biological Rationale

## Why TP53?

**TP53** encodes the tumour-suppressor protein p53, a major regulator of
cellular responses to genomic stress.

TP53 is involved in processes including:

- DNA-damage response;
- cell-cycle regulation;
- apoptosis;
- cellular senescence;
- genomic stability; and
- transcriptional regulation.

Recurrent TP53 mutations occur at specific positions in human cancers. These
mutation-associated hotspots provide a useful coordinate system for comparative
sequence analysis.

---

## Why elephants?

Elephants provide an important comparative system for evolutionary cancer
biology because they possess:

- large body size;
- long lifespan; and
- substantially more cells than humans.

Simple predictions based only on cell number and time at risk would suggest a
greater cancer burden in large, long-lived organisms.

The apparent lack of a straightforward scaling relationship is commonly
discussed in the context of **Peto's paradox**.

This project does **not** assume that TP53 sequence conservation alone explains
elephant cancer biology.

Instead, the computational question is narrower:

> **Do human TP53 mutation-associated positions exhibit conservation or
> sequence variation within elephant TP53-related sequences, and what
> comparative patterns can be identified for future investigation?**

---

# 🎯 Aim & Objectives

## Aim

To investigate the evolutionary conservation and sequence-level variation of
human TP53 mutation-associated hotspots across elephant TP53-related sequences
using comparative in-silico bioinformatics approaches.

## Objectives

- 🧬 Identify recurrent human TP53 mutation-associated hotspots.
- 🗂️ Retrieve and curate human and elephant TP53-related sequence resources.
- 🔎 Identify TP53-related elephant sequences using similarity-based analysis.
- 🧪 Perform protein-level sequence comparison and multiple sequence alignment.
- 🎯 Map human TP53 mutation-associated positions onto comparative sequences.
- 📊 Examine conservation and amino-acid variation at corresponding positions.
- 🐘 Compare canonical elephant TP53 with TP53-related retrogene sequence
  resources where applicable.
- 🌳 Establish an evolutionary context using phylogenetic analysis.
- 🤖 Explore sequence-derived features using computational modelling.
- ♻️ Maintain a transparent and reproducible research workflow.

---

# 🗂️ Data Sources

The project uses publicly available biological sequence resources.

## Human TP53

| Field | Information |
|---|---|
| Species | *Homo sapiens* |
| Gene | **TP53** |
| Database | UniProt |
| Accession | **P04637** |
| Sequence | Canonical protein |
| Local resource | `data/raw/human_tp53.fasta` |

Source:

**UniProt — TP53 P04637**

https://www.uniprot.org/uniprotkb/P04637/entry

---

## Elephant Genomic Resource

The principal elephant genomic resource documented for the project is:

| Field | Information |
|---|---|
| Species | *Elephas maximus* |
| Database | NCBI Assembly |
| Assembly | **GCA_024166365.1** |
| Assembly name | **mEleMax1 primary haplotype** |
| Role | Source for elephant TP53-related sequence identification |

Source:

**NCBI Assembly — GCA_024166365.1**

https://www.ncbi.nlm.nih.gov/assembly/13211691

---

## Large External Proteome

The original analysis used the corresponding elephant proteome, which was
approximately **41 MB**.

The complete proteome is not duplicated in this repository.

Instead, the repository preserves:

```text
Assembly accession
       ↓
Data provenance
       ↓
Derived TP53-related sequences
       ↓
Computational analysis
       ↓
Results
