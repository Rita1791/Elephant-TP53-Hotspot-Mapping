<a id="top"></a>

<div align="center">

🐘 Elephant TP53 Hotspot Mapping

The cancer-associated residues humans repeatedly mutate — traced through the elephant p53 family

Comparative genomics · Evolutionary cancer biology · Protein bioinformatics · Research software

Ritika Rajendra Rawat · Sermarani Nadar · Gursimran Kaur Uppal · 2026

<br>



<br>



<br>

<a href="figures/Figure7_composite_scores.png">
  <img src="figures/Figure7_composite_scores.png" alt="Human and elephant TP53 hotspot residue comparison" width="920">
</a>

<sub>👆 Click the hotspot comparison to open the full-resolution figure.</sub>

</div>

<div align="center">

🧭 Choose Your Route

🌍 ECCB 2026

🎯 Key Results

🧪 Research Journey

Follow the research progression

Explore the quantitative evidence

Trace how the analysis was built

🖼️ Figure Gallery

📂 Open Files

🤝 Connect

View full-resolution evidence

Jump directly into data and code

Contact the lead researcher

</div>

<a id="eccb-2026"></a>

🌍 From a Published Elephant Study to ECCB 2026

[!IMPORTANT]This focused human–elephant project is the published foundation of a broader mammalian TP53 research programme selected for poster presentation at ECCB 2026, the 25th European Conference on Computational Biology, in Geneva, Switzerland, from 31 August to 4 September 2026. Ms. Ritika Rajendra Rawat is the first and presenting author of the expanded study.

<div align="center">



</div>

Conference detail

Information

🏛️ Event

25th European Conference on Computational Biology — ECCB 2026

📍 Location

Geneva, Switzerland

📅 Dates

31 August–4 September 2026

🖼️ Format

Poster presentation

👩‍🔬 Presenter

Ms. Ritika Rajendra Rawat · first and presenting author

🧬 Expanded research

Evolutionary Conservation and Functional Constraint of TP53 Mutation Hotspots Across Mammalian Species

[!NOTE]The two repositories are connected stages of one research programme. They are not independent replications and not interchangeable datasets: this repository is human–elephant focused; the ECCB study uses an expanded 56-sequence mammalian dataset.

<p align="right"><a href="#top">⬆ Back to top</a></p>

⚡ The Study in 30 Seconds

This project asks one focused question:

How strongly are six recurrent human TP53 cancer hotspots preserved across curated elephant TP53 and TP53-related protein sequences?

The committed exploratory feature table contains 16 sequence records: one human TP53 reference and fifteen non-reference records. The hotspot table reports 90.91–100% exact-residue identity at R175, G245, R248, R249, R273 and R282. Protein-wide identity to human TP53 is much more variable, ranging from 52.42% to 83.33% among non-reference records.

That contrast is the scientific hook: strong local preservation at cancer-associated residues can coexist with substantial protein-wide divergence.

<table>
  <tr>
    <td align="center" width="25%"><h2>16</h2><sub>records in the committed<br>feature dataset</sub></td>
    <td align="center" width="25%"><h2>6</h2><sub>human TP53<br>hotspots evaluated</sub></td>
    <td align="center" width="25%"><h2>95.45%</h2><sub>mean hotspot identity<br>in the committed table</sub></td>
    <td align="center" width="25%"><h2>52.42–83.33%</h2><sub>non-reference identity<br>to human TP53</sub></td>
  </tr>
</table>

[!CAUTION]These are sequence-level, hypothesis-generating observations. They do not prove that TP53 alone explains elephant cancer resistance, resolve Peto’s paradox, or establish clinical utility.

<a id="key-results"></a>

🎯 Key Results — Explore the Evidence

<details open>
<summary><strong>🧬 Result 1 — All six selected hotspots show high exact-residue identity</strong></summary>

<br>

Human hotspot

Reference residue

Committed identity

Evidence state

R175

R

90.91%

🟢 Highly preserved

G245

G

100.00%

🟣 Fully preserved in the comparison

R248

R

90.91%

🟢 Highly preserved

R249

R

100.00%

🟣 Fully preserved in the comparison

R273

R

100.00%

🟣 Fully preserved in the comparison

R282

R

90.91%

🟢 Highly preserved

🔗 Open the evidence: hotspot CSV · focused hotspot figure · results summary

</details>

<details>
<summary><strong>📊 Result 2 — Local hotspot identity is stronger than global sequence similarity</strong></summary>

<br>

Mean identity across the six hotspot rows is 95.45%. Across the fifteen non-reference records in the committed feature table, mean protein-wide identity to human TP53 is 60.90%, with a median of 58.78%.

This is not a contradiction. Hotspot identity asks what happens at six human coordinates; global identity measures similarity across the aligned protein.

🔗 Open the evidence: feature table · identity plot · methodology

</details>

<details>
<summary><strong>🏆 Result 3 — Two non-reference records lead the similarity ranking</strong></summary>

<br>

XP_049714738.1 and XP_003416950.3 each show 83.33% identity to the human TP53 reference in the committed exploratory table. This ranking is descriptive; it does not by itself prove orthology, biochemical equivalence or functional activity.

🔗 Open the evidence: top-similarity table · complete feature table

</details>

<details>
<summary><strong>🧩 Result 4 — The repository contains distinct evidence layers</strong></summary>

<br>

Evidence layer

Scope

Correct use

🖼️ Focused project figures

Selected canonical and TP53-related examples

Visual residue correspondence in the displayed comparison

📊 Historical results/ML/ tables

16 feature records; 15 non-reference records

Quantitative exploratory results already committed

🐍 Current analysis script

Explicit pairwise coordinate mapping and newer outputs

Reusable implementation that still needs reconciliation with historical artifacts

The focused figure and the broader quantitative table are related, but they do not use the same denominator. Every numerical claim should therefore link to the exact artifact that supports it.

🔗 Read the evidence contract: reproducibility · artifact index

</details>

<details>
<summary><strong>🐘 Result 5 — The analysis became an explorable research prototype</strong></summary>

<br>

EleProtect v2.0 turns the sequence-analysis workflow into a Streamlit interface for TP53-oriented preprocessing, hotspot mapping, feature extraction, exploratory ranking, visualisation and CSV export.

It is research software—not a diagnostic, prognostic or treatment-selection system.

🔗 Explore EleProtect: application · app guide · source code

</details>

<p align="right"><a href="#top">⬆ Back to top</a></p>

<a id="research-journey"></a>

🧪 Research Journey

flowchart TB
    A["🧑 Human TP53 · P04637"] --> B["🐘 Curated elephant TP53-related proteins"]
    B --> C["🧹 Validate and preprocess sequences"]
    C --> D["🧬 Align proteins to the human reference"]
    D --> E["🎯 Map six cancer-associated hotspots"]
    E --> F["📊 Measure local and global similarity"]
    F --> G["🌳 Add clustering and phylogenetic context"]
    G --> H["🐘 Explore results through EleProtect"]
    H --> I["⚖️ Interpret without causal or clinical overclaiming"]

🔎 Follow Each Stage into the Repository

Stage

What happens

Open the file

01 · 🧬

Load the analysis-ready human–elephant protein dataset

TP53_clean.fasta

02 · ✅

Validate sequence symbols, identifiers and human reference

analysis script

03 · 🧭

Align comparative sequences to human TP53

methodology

04 · 🎯

Map R175, G245, R248, R249, R273 and R282

hotspot table

05 · 📏

Calculate sequence features and global identity

feature table

06 · 🤖

Explore descriptive feature-space clusters

clustered features

07 · 🌳

Place sequences in phylogenetic context

Newick tree

08 · 🐘

Expose the workflow through a research interface

EleProtect app

09 · ⚖️

Separate observed evidence from causal claims

limitations

<details>
<summary><strong>🔬 Expand the complete analytical design</strong></summary>

Curate human and elephant TP53-related protein resources with source traceability.

Clean sequences and verify valid amino-acid symbols.

Anchor the comparison to canonical human TP53, UniProt P04637.

Perform explicit global pairwise alignment.

Map one-based human coordinates through each alignment.

Classify each focal position as identical, substituted, gapped or unmapped.

Exclude the human reference from non-reference similarity summaries.

Calculate sequence features and global identity to human TP53.

Use clustering as a descriptive feature-space view—not a biological classifier.

Preserve machine-readable outputs, figures and interpretation boundaries.

📚 Full documentation: methodology · provenance · reproducibility

</details>

<p align="right"><a href="#top">⬆ Back to top</a></p>

<a id="figure-gallery"></a>

🖼️ Interactive Figure Gallery

Click any figure to open the original high-resolution repository image.

<table>
  <tr>
    <td width="50%" align="center">
      <a href="figures/Figure7_composite_scores.png">
        <img src="figures/Figure7_composite_scores.png" alt="Human and elephant TP53 hotspot residue comparison" width="100%">
      </a><br>
      <strong>01 · 🎯 Focused hotspot comparison</strong><br>
      <sub>Human, African elephant, Asian elephant and retrogene residue correspondence.</sub>
    </td>
    <td width="50%" align="center">
      <a href="results/ML/identity_barplot.png">
        <img src="results/ML/identity_barplot.png" alt="Identity of TP53-related records to human TP53" width="100%">
      </a><br>
      <strong>02 · 📊 Global sequence identity</strong><br>
      <sub>Protein-wide similarity tells a broader, more variable story.</sub>
    </td>
  </tr>
  <tr>
    <td width="50%" align="center">
      <a href="figures/Figure4_MSA.png">
        <img src="figures/Figure4_MSA.png" alt="Human and elephant TP53 multiple-sequence alignment" width="100%">
      </a><br>
      <strong>03 · 🧬 Alignment landscape</strong><br>
      <sub>Sequence-level context behind the hotspot mapping.</sub>
    </td>
    <td width="50%" align="center">
      <a href="figures/Figure6_phylogenetic_tree.png">
        <img src="figures/Figure6_phylogenetic_tree.png" alt="Phylogenetic tree of analysed TP53-related records" width="100%">
      </a><br>
      <strong>04 · 🌳 Phylogenetic context</strong><br>
      <sub>Relationships among the analysed sequence records.</sub>
    </td>
  </tr>
  <tr>
    <td width="50%" align="center">
      <a href="figures/Figure2_preprocessing_workflow.png">
        <img src="figures/Figure2_preprocessing_workflow.png" alt="TP53 sequence preprocessing workflow" width="48%">
      </a><br>
      <strong>05 · 🧹 Sequence preprocessing</strong><br>
      <sub>From retrieved protein records to analysis-ready inputs.</sub>
    </td>
    <td width="50%" align="center">
      <a href="figures/Figure3_pipeline.png">
        <img src="figures/Figure3_pipeline.png" alt="EleProtect analysis pipeline" width="48%">
      </a><br>
      <strong>06 · 🐘 EleProtect workflow</strong><br>
      <sub>The application-facing analysis and exploration path.</sub>
    </td>
  </tr>
</table>

<div align="center">



</div>

<p align="right"><a href="#top">⬆ Back to top</a></p>

🐘 EleProtect v2.0

<table>
  <tr>
    <td width="24%" align="center">
      <a href="EleProtect_App/">
        <img src="EleProtect_App/EleProtect.png" alt="EleProtect v2.0 research software" width="165">
      </a>
    </td>
    <td width="76%">
      <h3>Move from static outputs to interactive exploration</h3>
      EleProtect is the Streamlit research interface associated with this project. It supports TP53-oriented sequence processing, hotspot mapping, feature extraction, exploratory ranking, visualisation and CSV export.<br><br>
      <strong>Boundary:</strong> research prototype only; not validated for clinical decision-making.
    </td>
  </tr>
</table>

<div align="center">



</div>

<a id="open-the-research"></a>

📂 Open the Research

Pick What You Need

I want to…

Start here

Direct route

⚡ Understand the project quickly

Reviewer summary

docs/reviewer_summary.md

🎯 Inspect the quantitative findings

Results summary

results/summary.md

🧪 Examine the complete method

Methodology

docs/methodology.md

🧾 Verify sequence sources

Provenance

docs/provenance.md

🧬 Inspect analysis-ready sequences

Processed FASTA

data/processed/TP53_clean.fasta

💻 Audit the implementation

Main Python workflow

scripts/TP53_Comparative_Analysis.py

📓 Explore the core analysis

Jupyter notebook

TP53_Comparative_Analysis.ipynb

🤖 Inspect feature prioritisation

Python workflow

EleProtect_Feature_Prioritization.py

🌳 Reuse the phylogeny

Newick tree

results/phylogeny/TP53_tree.nwk

🖼️ View the project poster

Poster PDF

docs/Srijna_Poster.pdf

📚 Read the published chapter

Publication PDF

manuscript/Published_Research.pdf

✍️ Cite the project

Citation metadata

CITATION.cff

⚠️ Check limitations first

Interpretation boundaries

docs/interpretation_and_limitations.md

Repository Map

Elephant-TP53-Hotspot-Mapping/
├── 🧬 data/             Curated and analysis-ready sequence resources
├── 🐍 scripts/          Reusable comparative-analysis implementation
├── 📓 notebooks/        Interactive and exploratory analyses
├── 📊 results/          Alignments, feature tables, trees and summaries
├── 🖼️ figures/          Research diagrams and visual evidence
├── 🐘 EleProtect_App/   Streamlit research prototype
├── 📚 docs/             Methods, provenance, review notes and limitations
├── 📘 manuscript/       Published research and associated material
├── ✅ tests/            Unit and result-contract checks
└── ✍️ CITATION.cff      Machine-readable citation metadata

<details>
<summary><strong>💻 Clone and run the current analysis</strong></summary>

git clone https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping.git
cd Elephant-TP53-Hotspot-Mapping

python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

python scripts/TP53_Comparative_Analysis.py
pytest -q

[!CAUTION]The current script writes a newer output contract than the historical exploratory tables under results/ML/. Until those layers are fully reconciled, cite the exact artifact used for each result.

</details>

<p align="right"><a href="#top">⬆ Back to top</a></p>

🔍 Scientific Boundaries

✅ Supported by this repository

❌ Not established by this repository

Mapping six human TP53 hotspots onto curated elephant-related sequences

That TP53 alone explains elephant cancer resistance

High exact-residue identity at the selected positions

A complete molecular solution to Peto’s paradox

Protein-wide identity and sequence-feature comparisons

Functional activity of every TP53-related sequence

Descriptive clustering and phylogenetic context

Clinical risk, diagnosis, prognosis or treatment response

A traceable foundation for broader mammalian research

Experimental validation of the proposed biological interpretation

<details>
<summary><strong>⚠️ Read the limitations and reproducibility status</strong></summary>

Exact amino-acid identity is not a codon-aware evolutionary model. Results depend on sequence selection, annotation quality and alignment decisions. Protein similarity does not demonstrate biochemical equivalence. Feature clusters are exploratory and are not evolutionary lineages. BLOSUM62 values are substitution scores—not pathogenicity estimates.

The current implementation and historical committed tables also use different output contracts. The repository supports evidence inspection and modular analysis, but it should not be advertised as a fully reconciled one-command rebuild until those layers are aligned.

🔗 Read: full limitations · reproducibility · project relationship

</details>

<details>
<summary><strong>🧬 Understand the elephant and mammalian projects</strong></summary>

Research stage

Primary question

Correct use

🐘 Human–elephant study

How do human hotspots map onto elephant TP53-related sequences?

Published foundation, pairwise mapping and EleProtect

🌍 Expanded mammalian study

Are recurrent human hotspots unusually constrained across mammals?

56-sequence conservation, cancer recurrence and phylogenetic sensitivity

The mammalian repository extends the biological question. It does not retroactively change the data or results of this focused elephant project.

</details>

📖 Publication & Citation

Rawat, R. R., Nadar, S., & Uppal, G. K. (2026).Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants: A Responsible Bioinformatics Innovation Contributing to Cancer Research.In Ideathon on Vikshit Bharat: Ideas, Innovation and Impact, Chapter 11, pp. 123–138. https://doi.org/10.25215/9141002199

<div align="center">



</div>

<details>
<summary><strong>📋 Copy the BibTeX citation</strong></summary>

@incollection{rawat2026elephanttp53,
  author    = {Rawat, Ritika Rajendra and Nadar, Sermarani and Uppal, Gursimran Kaur},
  title     = {Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants: A Responsible Bioinformatics Innovation Contributing to Cancer Research},
  booktitle = {Ideathon on Vikshit Bharat: Ideas, Innovation and Impact},
  chapter   = {11},
  pages     = {123--138},
  year      = {2026},
  publisher = {RED'SHINE Publication},
  doi       = {10.25215/9141002199}
}

</details>

👩‍🔬 Research Team

Researcher

Contribution

Ritika Rajendra Rawat

Study conception, comparative workflow, sequence analysis, research software, interpretation, visualisation and scientific communication

Sermarani Nadar

Scientific discussion, interpretation, critical review and co-authorship of the associated publication

Gursimran Kaur Uppal

Scientific discussion, interpretation, critical review and co-authorship of the associated publication

<a id="connect-with-the-researcher"></a>

🤝 Connect with the Researcher

<div align="center">

Ritika Rajendra Rawat

Bioinformatics researcher working across comparative genomics, evolutionary cancer biology, computational biology and evidence-aware genomic analysis.

Interested in research collaboration, ECCB discussion, TP53 biology, PhD-related academic communication or reproducible bioinformatics? Choose the route that matches your purpose:

<br>



<br>

Purpose

Direct action

🤝 Research collaboration

Send a collaboration email

💼 Professional networking

Connect on LinkedIn

💻 Code and project portfolio

Visit GitHub profile

🔬 Publications and research updates

Follow on ResearchGate

🐛 Technical question about this repository

Open a GitHub issue

</div>

📜 License

Code and repository-authored material are released under the MIT License. Third-party biological data and external resources retain their original terms of use.

<div align="center">

🐘 The elephant is the biological question.

🧬 The conserved residue is the evidence trail.

🌍 ECCB 2026 · 🎯 Results · 🧪 Journey · 🖼️ Figures · 📂 Files · 🤝 Connect · ⬆ Back to top

<sub>Comparative-genomics research · Hypothesis-generating · Not validated for clinical use</sub>

</div>
