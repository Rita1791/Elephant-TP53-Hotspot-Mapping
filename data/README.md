# 🧬 Research Data

This directory contains the biological sequence resources, processed datasets,
and derived reference tables used throughout the **Elephant TP53 Hotspot Mapping**
research project.

The data organization follows a simple computational-research principle:

```text
PUBLIC BIOLOGICAL RESOURCES
          │
          ▼
       raw/
          │
          ▼
    preprocessing
          │
          ▼
    processed/
          │
          ▼
 comparative analysis
          │
          ▼
     results/

data/
│
├── raw/
│   ├── Blast/
│   ├── ncbi_dataset/
│   ├── african_tp53_hits.fasta
│   ├── asian_tp53_hits.fasta
│   ├── elephant_assembly_data_report.jsonl
│   ├── elephant_tp53_RTG9.fasta
│   ├── elephant_tp53_RTG9_translate.fasta
│   ├── elephant_tp53_canonical.fasta
│   ├── human_tp53.fasta
│   └── provenance.txt
│
├── processed/
│   ├── TP53_all_sequences.fasta
│   ├── TP53_clean.fasta
│   ├── human_elephant_tp53_pair.fasta
│   └── human_elephant_tp53_retrogene_comparison.fasta
│
└── Database/
    ├── README.md
    ├── tp53_elephant_database.csv
    └── tp53_elephant_database.json

                    PUBLIC DATABASES
                          │
                          ▼
                 ┌─────────────────┐
                 │     raw/        │
                 │                 │
                 │ NCBI / BLAST /  │
                 │ sequence files  │
                 └────────┬────────┘
                          │
                          ▼
                  sequence curation
                          │
                          ▼
                 ┌─────────────────┐
                 │   processed/    │
                 │                 │
                 │ cleaned / paired│
                 │ comparative     │
                 │ sequences       │
                 └────────┬────────┘
                          │
                          ▼
                 comparative analysis
                          │
             ┌────────────┼────────────┐
             ▼            ▼            ▼
           MSA        Hotspot       Similarity
                      Mapping        Analysis
             │            │            │
             └────────────┼────────────┘
                          ▼
                       results/
                          │
                          ▼
                    Interpretation
                          │
                          ▼
                  Research Outputs

