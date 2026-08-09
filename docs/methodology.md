# 🧬 Methodology

<p align="center">

<b>Comparative in-silico framework for TP53 sequence analysis between humans and elephants</b>

</p>

<p align="center">

🧬 Sequence Analysis &nbsp; • &nbsp;
🔎 Similarity Analysis &nbsp; • &nbsp;
🧪 Multiple Sequence Alignment &nbsp; • &nbsp;
📍 Hotspot Mapping &nbsp; • &nbsp;
🌳 Phylogenetics &nbsp; • &nbsp;
🤖 Feature Analysis

</p>

---

## 🎯 1. Study Objective

The **Elephant TP53 Hotspot Mapping** project uses a comparative
bioinformatics framework to investigate TP53-related sequence characteristics
between humans and elephants.

The analysis focuses on:

- identification of TP53-related elephant sequences;
- comparison of human and elephant TP53 sequences;
- sequence conservation and divergence;
- mapping of human TP53 mutation-associated positions;
- evolutionary relationships among analysed sequences;
- extraction of sequence-derived features; and
- exploratory computational prioritization.

The overall objective is to identify **comparative sequence patterns that may
provide candidates for further biological investigation**.

The analysis is computational and does not independently establish experimental
mechanisms.

---

# 🧭 2. Overall Computational Workflow

The project follows the workflow:

```text
                    REFERENCE DATA
                         │
                         ▼
              ┌─────────────────────┐
              │ Sequence Retrieval  │
              └──────────┬──────────┘
                         │
                         ▼
              ┌─────────────────────┐
              │ Data Preprocessing  │
              └──────────┬──────────┘
                         │
                         ▼
              ┌─────────────────────┐
              │ Similarity Analysis │
              └──────────┬──────────┘
                         │
                         ▼
              ┌─────────────────────┐
              │ Sequence Selection  │
              └──────────┬──────────┘
                         │
                         ▼
              ┌─────────────────────┐
              │ Multiple Sequence   │
              │ Alignment           │
              └──────────┬──────────┘
                         │
              ┌──────────┴──────────┐
              ▼                     ▼
     ┌─────────────────┐   ┌─────────────────┐
     │ Hotspot Mapping │   │ Phylogenetic    │
     │                 │   │ Analysis        │
     └────────┬────────┘   └────────┬────────┘
              │                     │
              └──────────┬──────────┘
                         ▼
              ┌─────────────────────┐
              │ Feature Extraction  │
              └──────────┬──────────┘
                         │
                         ▼
              ┌─────────────────────┐
              │ Exploratory         │
              │ Computational       │
              │ Modelling           │
              └──────────┬──────────┘
                         │
                         ▼
              ┌─────────────────────┐
              │ Visualization &     │
              │ Interpretation      │
              └─────────────────────┘
