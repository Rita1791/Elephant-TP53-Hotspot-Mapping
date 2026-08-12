<a id="top"></a>

<div align="center">

<img src="EleProtect_App/EleProtect.png" alt="EleProtect — Elephant TP53 research emblem" width="145">

🐘 Elephant TP53 Hotspot Mapping 🧬

Tracing human cancer-mutation hotspots through the elephant p53 sequence landscape

Comparative genomics · Evolutionary cancer biology · Protein sequence analysis · Reproducible research

Ritika Rajendra Rawat · Sermarani Nadar · Gursimran Kaur Uppal · 2026

<br>



<br>



<br>

Six human TP53 hotspots. Fifteen non-reference sequences. One evolutionary question.

What remains conserved when the surrounding protein landscape changes?

</div>

<div align="center">

🧭 Choose Your Route

⚡ 60-Second Tour

🎯 Key Results

🧪 Research Workflow

Understand the project

Inspect quantitative evidence

Follow data to interpretation

🖼️ Figure Gallery

📂 Open the Research

🤝 Connect

Explore visual outputs

Jump into files and code

Contact the lead researcher

</div>

<a id="eccb-2026"></a>

🌍 From Elephant TP53 to ECCB 2026

[!IMPORTANT]This elephant-focused project is the published foundation of a broader mammalian TP53 research programme selected for poster presentation at ECCB 2026, the 25th European Conference on Computational Biology. Ms. Ritika Rajendra Rawat is the first and presenting author.

<table>
  <tr>
    <td width="50%">
      <strong>🖼️ Accepted poster</strong><br><br>
      <em>Evolutionary Conservation and Functional Constraints of TP53 Mutation Hotspots Across Mammalian Species</em>
    </td>
    <td width="50%">
      <strong>📍 Geneva, Switzerland</strong><br><br>
      📅 31 August–4 September 2026<br>
      🎫 Poster C-G.32 · Session 3<br>
      👩‍🔬 Presenter: Ritika Rajendra Rawat
    </td>
  </tr>
</table>

<div align="center">



</div>

Research-lineage boundary: ECCB 2026 presents the expanded mammalian analysis. This repository contains the earlier, focused human–elephant study and EleProtect prototype. The two datasets should not be described as independent replications.

<p align="right"><a href="#top">⬆ Back to top</a></p>

<a id="overview"></a>

⚡ The Research in 60 Seconds

Human cancers repeatedly alter a small set of functionally important TP53 residues. This project maps six of those human hotspots onto curated elephant TP53 and TP53-related protein sequences, then asks whether local residue conservation persists despite broader sequence divergence.

<table>
  <tr>
    <td align="center" width="25%"><h2>16</h2><sub>sequence records in the<br>committed feature dataset</sub></td>
    <td align="center" width="25%"><h2>6</h2><sub>human TP53<br>hotspots evaluated</sub></td>
    <td align="center" width="25%"><h2>95.45%</h2><sub>mean hotspot identity in the<br>committed summary</sub></td>
    <td align="center" width="25%"><h2>52.42–83.33%</h2><sub>non-reference identity<br>to human TP53</sub></td>
  </tr>
</table>

The headline

All six selected hotspots were highly conserved in the committed comparison: G245, R249 and R273 reached 100%, while R175, R248 and R282 reached 90.91%.

The honest interpretation

This is sequence-level computational evidence consistent with evolutionary constraint. It does not prove that TP53 conservation causes elephant cancer resistance, explain Peto's paradox by itself, or establish clinical utility.

🔗 Inspect the source: results/summary.md · results/ML/TP53_hotspot_analysis.csv · docs/interpretation_and_limitations.md

<a id="question"></a>

🧠 Why Elephants? Why TP53?

Elephants combine a large body, long lifespan and enormous number of lifetime cell divisions—conditions that should increase opportunities for malignant transformation. Their cancer burden does not scale as simply as those factors predict. That tension is commonly discussed as Peto's paradox.

flowchart TD
    A["🐘 Large body + long lifespan"] --> B["Many lifetime cell divisions"]
    B --> C["Expected higher cancer burden"]
    C --> D["Peto's paradox"]
    D --> E["🧬 Compare TP53 sequence constraints"]

TP53 is not assumed to be the complete answer. It is used as a focused comparative system because p53 is central to DNA-damage responses, cell-cycle control, senescence and apoptosis—and because recurrent human cancer mutations cluster at defined residues.

<details>
<summary><strong>🎯 Open the biological question and study objective</strong></summary>

Research question

How do canonical human TP53 cancer-mutation hotspot residues map onto curated elephant TP53 and TP53-related protein sequences?

Aim

To investigate sequence-level conservation and variation at major human TP53 mutation-associated hotspots across elephant TP53-related sequences using comparative in-silico bioinformatics.

Six focal positions

R175 · G245 · R248 · R249 · R273 · R282

The human reference is canonical TP53 protein UniProt P04637, 393 amino acids long.

🔗 Open UniProt P04637 · Read the reviewer summary

</details>

<p align="right"><a href="#top">⬆ Back to top</a></p>

<a id="workflow"></a>

🧪 Research Workflow

flowchart TD
    A["🧑 Human TP53 · P04637"] --> B["🐘 Elephant TP53-related sequences"]
    B --> C["🧹 Curate + validate proteins"]
    C --> D["🧬 Global sequence alignment"]
    D --> E["🎯 Map six human hotspots"]
    E --> F{"Two analytical lenses"}
    F --> G["Local residue identity"]
    F --> H["Global sequence features"]
    G --> I["Conservation evidence"]
    H --> J["Similarity + clustering"]
    I --> K["⚖️ Responsible interpretation"]
    J --> K

Follow the evidence into the repository

Stage

What happens

Open the artifact

01 · 🧬

Load the curated human–elephant protein dataset

TP53_clean.fasta

02 · ✅

Validate protein symbols, IDs and the human reference

analysis script

03 · 🧭

Align each comparative sequence to human TP53

methodology

04 · 🎯

Map R175, G245, R248, R249, R273 and R282

hotspot table

05 · 📏

Measure broader identity and sequence features

feature table

06 · 🤖

Explore similarity-based prioritisation and clustering

clustered features

07 · 🌳

Place the analysed sequences in phylogenetic context

Newick tree

08 · ⚖️

Separate computational observations from causal claims

limitations

<details>
<summary><strong>🔬 Expand the full analytical design</strong></summary>

Retrieve human and elephant TP53-related sequence resources.

Preserve accession and source information under data/metadata/ and docs/provenance.md.

Clean and validate protein sequences.

Identify the human reference using UniProt accession P04637.

Perform explicit global pairwise alignment.

Map one-based human TP53 coordinates through each alignment.

Classify each mapped position as identical, substituted, gapped or unmapped.

Exclude the human reference from comparative denominators.

Calculate sequence features and global identity to human TP53.

Explore feature-space clustering as a descriptive—not biological—classifier.

Preserve machine-readable tables, figures and interpretation boundaries.

📚 Full documentation: methodology · provenance · reproducibility

</details>

<p align="right"><a href="#top">⬆ Back to top</a></p>

<a id="key-results"></a>

🎯 Key Results — Explore the Evidence

1 · Hotspot identity

Human hotspot

Reference residue

Committed result

Evidence state

R175

R

90.91%

🟢 Highly conserved

G245

G

100.00%

🟣 Fully conserved

R248

R

90.91%

🟢 Highly conserved

R249

R

100.00%

🟣 Fully conserved

R273

R

100.00%

🟣 Fully conserved

R282

R

90.91%

🟢 Highly conserved

🔗 Open the machine-readable hotspot CSV

2 · Global similarity tells a different story

Local hotspot preservation coexists with much broader protein-level divergence.

<table>
  <tr>
    <td align="center" width="25%"><h3>83.33%</h3><sub>highest non-reference<br>identity</sub></td>
    <td align="center" width="25%"><h3>60.90%</h3><sub>mean non-reference<br>identity</sub></td>
    <td align="center" width="25%"><h3>58.78%</h3><sub>median non-reference<br>identity</sub></td>
    <td align="center" width="25%"><h3>3</h3><sub>exploratory feature<br>clusters</sub></td>
  </tr>
</table>

Two records—XP_049714738.1 and XP_003416950.3—show the highest reported non-reference identity to human TP53 at 83.33%.

🔗 Open the similarity ranking · Open clustered features

3 · Understand the two evidence layers

Evidence layer

Scope

Correct interpretation

🖼️ Focused project figures

Selected canonical and TP53-related examples

Visual residue correspondence in the displayed comparison

📊 Committed results/ML/ tables

16 feature records; 15 non-reference records

Historical quantitative exploratory summary

🐍 Current analysis script

Explicit pairwise mapping with exact status fields

Reproducible implementation that generates the modern output contract when run

[!NOTE]Do not treat the focused figure and the broader quantitative denominator as the same comparison. The README keeps both because they document different stages of the project.

<details>
<summary><strong>🧩 Why local conservation and global identity can differ</strong></summary>

Global sequence identity measures similarity across an aligned protein. Hotspot identity asks a narrower question at six human reference coordinates. A sequence can therefore diverge substantially overall while retaining residues at locally constrained positions.

That pattern supports further evolutionary and structural investigation; it does not prove shared protein function.

</details>

<a id="figures"></a>

🖼️ Interactive Figure Gallery

Click any figure to open the original, full-resolution repository image.

<a href="results/ML/identity_barplot.png">
  <img src="results/ML/identity_barplot.png" alt="Identity of elephant TP53-related sequences to human TP53" width="100%">
</a>

<p align="center"><strong>📊 Global sequence identity to human TP53</strong><br><sub>Human TP53 is the 100% reference; non-reference records span 52.42–83.33% in the committed feature dataset.</sub></p>

<table>
  <tr>
    <td width="50%" align="center">
      <a href="figures/Figure7_composite_scores.png">
        <img src="figures/Figure7_composite_scores.png" alt="Focused TP53 hotspot mapping table" width="100%">
      </a><br>
      <strong>01 · 🎯 Focused hotspot map</strong><br>
      <sub>Residue correspondence across selected human and elephant sequences.</sub>
    </td>
    <td width="50%" align="center">
      <a href="figures/Figure4_MSA.png">
        <img src="figures/Figure4_MSA.png" alt="Multiple sequence alignment of TP53-related proteins" width="100%">
      </a><br>
      <strong>02 · 🧬 Multiple sequence alignment</strong><br>
      <sub>The alignment context behind residue-level comparison.</sub>
    </td>
  </tr>
  <tr>
    <td width="50%" align="center">
      <a href="figures/Figure6_phylogenetic_tree.png">
        <img src="figures/Figure6_phylogenetic_tree.png" alt="Phylogenetic tree of analysed TP53-related sequences" width="100%">
      </a><br>
      <strong>03 · 🌳 Phylogenetic context</strong><br>
      <sub>Sequence relationships in the analysed dataset.</sub>
    </td>
    <td width="50%" align="center">
      <a href="figures/Figure3_pipeline.png">
        <img src="figures/Figure3_pipeline.png" alt="Elephant TP53 computational analysis pipeline" width="50%">
      </a><br>
      <strong>04 · ⚙️ Analysis pipeline</strong><br>
      <sub>From sequence resources to comparative interpretation.</sub>
    </td>
  </tr>
</table>

<div align="center">



</div>

<p align="right"><a href="#top">⬆ Back to top</a></p>

<a id="eleprotect"></a>

💻 EleProtect v2.0 — Research Software

<table>
  <tr>
    <td width="25%" align="center">
      <a href="EleProtect_App/">
        <img src="EleProtect_App/EleProtect.png" alt="EleProtect v2.0" width="170">
      </a>
    </td>
    <td width="75%">
      <h3>Sequence analysis made explorable</h3>
      EleProtect is the Streamlit research interface associated with this project. It accepts DNA or protein sequence input, performs TP53-oriented processing and feature extraction, exposes exploratory ML ranking, and supports CSV export.<br><br>
      <strong>Research prototype only:</strong> it is not a diagnostic, prognostic or therapeutic decision-support system.
    </td>
  </tr>
</table>

<div align="center">



</div>

cd EleProtect_App
python -m pip install -r requirements.txt
streamlit run app.py

<a id="research-explorer"></a>

📂 Open the Research

Pick what you need

I want to…

Best starting point

Direct route

⚡ Understand the project quickly

Reviewer summary

docs/reviewer_summary.md

🎯 Inspect the quantitative findings

Results summary

results/summary.md

🧪 Examine the method

Methodology

docs/methodology.md

🧾 Verify sequence sources

Provenance

docs/provenance.md

🧬 Inspect the analysis-ready sequences

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

🖼️ View the student poster

Srijna poster

docs/Srijna_Poster.pdf

📚 Read the published chapter

Publication PDF

manuscript/Published_Research.pdf

⚠️ Check scientific boundaries

Interpretation & limitations

docs/interpretation_and_limitations.md

✍️ Cite the project

Citation metadata

CITATION.cff

Repository map

Elephant-TP53-Hotspot-Mapping/
├── 🧬 data/             Raw, processed and metadata-backed sequence resources
├── 🐍 scripts/          Reusable comparative-analysis implementation
├── 📓 notebooks/        Interactive and exploratory computational work
├── 📊 results/          Alignments, feature tables, trees and summaries
├── 🖼️ figures/          Research diagrams and visual outputs
├── 🐘 EleProtect_App/   Streamlit research prototype
├── 📚 docs/             Methods, provenance, review guide and limitations
├── 📘 manuscript/       Thesis and associated published research
├── ✅ tests/            Unit and result-contract checks
└── ✍️ CITATION.cff      Machine-readable citation metadata

<a id="reproduce"></a>

♻️ Reproduce the Current Analysis

git clone https://github.com/Rita1791/Elephant-TP53-Hotspot-Mapping.git
cd Elephant-TP53-Hotspot-Mapping

python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

python scripts/TP53_Comparative_Analysis.py
pytest -q

<details>
<summary><strong>🪟 Windows activation command</strong></summary>

.venv\Scripts\Activate.ps1

</details>

<details>
<summary><strong>📦 Outputs generated by the current script</strong></summary>

results/tp53_sequence_features.csv
results/tp53_hotspot_mapping.csv
results/tp53_hotspot_identity_summary.csv
results/tp53_comparative_features.csv
results/tp53_comparative_features_clustered.csv
results/tp53_excluded_sequences.csv
results/tp53_summary.json

figures/tp53_identity_to_human.png
figures/tp53_hotspot_conservation.png
figures/tp53_feature_clustering.png

</details>

[!CAUTION]The current script writes a newer output contract than the historical exploratory tables under results/ML/. Until the regenerated outputs are committed and reconciled with the older artifacts, cite the exact file used for each numerical claim.

<a id="boundaries"></a>

⚖️ What the Evidence Supports

✅ Supported by this repository

❌ Not established by this repository

Mapping six human TP53 hotspots onto curated elephant-related sequences

That TP53 alone explains elephant cancer resistance

High exact residue identity at the selected positions

A complete mechanistic solution to Peto's paradox

Global identity and sequence-feature comparisons

Functional activity of every TP53-related sequence

Descriptive phylogenetic and clustering context

Clinical risk, diagnosis, prognosis or treatment response

A traceable foundation for follow-up research

Experimental validation of molecular function

<details>
<summary><strong>⚠️ Read the scientific and computational limitations</strong></summary>

Results depend on sequence selection, annotation quality and alignment parameters.

Exact amino-acid identity is not a formal phylogenetic conservation score.

Sequence similarity does not demonstrate shared biochemical activity.

The feature-based clusters are exploratory and are not evolutionary lineages.

The BLOSUM62 values are substitution scores, not pathogenicity or clinical-impact estimates.

The elephant-focused and expanded mammalian projects answer related but different questions.

🔗 Full limitations · reproducibility record · companion-repository relationship

</details>

<a id="publication"></a>

📖 Publication & Citation

Rawat, R. R., Nadar, S., & Uppal, G. K. (2026).Comparative In-Silico Mapping of TP53 Mutation Hotspots in Elephants: A Responsible Bioinformatics Innovation Contributing to Cancer Research.In Ideathon on Vikshit Bharat: Ideas, Innovation and Impact, Chapter 11, pp. 123–138.https://doi.org/10.25215/9141002199

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

👥 Research Team

Researcher

Contribution

Ritika Rajendra Rawat

Primary researcher; comparative workflow, sequence analysis, research development, software and scientific communication

Sermarani Nadar

Co-author of the associated published human–elephant TP53 research

Gursimran Kaur Uppal

Co-author of the associated published human–elephant TP53 research

<a id="connect"></a>

🤝 Connect with the Researcher

<div align="center">

👩‍🔬 Ritika Rajendra Rawat

MSc Bioinformatics · Researcher working across comparative genomics, evolutionary cancer biology, computational biology and evidence-aware genomic analysis.

For research collaboration, TP53 discussion, PhD-related academic communication or technical questions, choose the direct route that fits your purpose:

<br>



<br>

Purpose

Direct action

🤝 Research collaboration

Send a collaboration email

💼 Professional networking

Connect on LinkedIn

💻 Code and project portfolio

Visit the GitHub profile

🔬 Publications and research updates

Follow on ResearchGate

🐛 Technical question about this project

Open a GitHub issue

</div>

📜 License

Code and repository-authored materials are released under the MIT License. External biological data, software and third-party resources retain their original terms.

<div align="center">

🐘 The elephant is the biological question.

🧬 The conserved residue is the evidence trail.

⚡ Overview · 🎯 Results · 🧪 Workflow · 🖼️ Figures · 📂 Files · 🤝 Connect · ⬆ Back to top

<sub>Comparative bioinformatics research · Hypothesis-generating · Not validated for clinical use</sub>

</div>
