{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# 🐘 EleProtect — Machine Learning Analysis\n",
    "\n",
    "## Exploratory Computational Framework for Elephant TP53 Sequence Features\n",
    "\n",
    "**Research Project:** Elephant TP53 Hotspot Mapping  \n",
    "**Researcher:** Ritika Rajendra Rawat  \n",
    "**Degree:** MSc Bioinformatics  \n",
    "**Institution:** University of Mumbai\n",
    "\n",
    "---\n",
    "\n",
    "## 🎯 Objective\n",
    "\n",
    "This notebook provides an exploratory machine-learning analysis of sequence-derived features associated with the EleProtect research framework.\n",
    "\n",
    "The objective is to characterize protein-sequence features computationally and establish a transparent foundation for future statistical or machine-learning studies.\n",
    "\n",
    "The analysis includes:\n",
    "\n",
    "- sequence quality assessment;\n",
    "- amino-acid composition;\n",
    "- physicochemical sequence descriptors;\n",
    "- feature normalization;\n",
    "- exploratory dimensionality reduction;\n",
    "- unsupervised clustering;\n",
    "- feature correlation analysis;\n",
    "- exploratory feature ranking;\n",
    "- reproducible result export.\n",
    "\n",
    "> ⚠️ This notebook is an exploratory computational analysis. It is not a clinical model, cancer-risk predictor, cancer-resistance predictor, pathogenicity classifier, or diagnostic system."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# 🧬 Scientific Rationale\n",
    "\n",
    "Protein sequences contain measurable compositional and physicochemical properties that can be represented numerically for computational analysis.\n",
    "\n",
    "Characterizing these properties can help identify patterns, similarities, and differences among TP53-related sequences and can provide candidate features for future research.\n",
    "\n",
    "The present workflow is deliberately transparent and exploratory. Feature-level differences should be treated as descriptive observations rather than direct evidence of biological function.\n",
    "\n",
    "Any future predictive model would require appropriate biological labels, sufficiently large datasets, independent validation, baseline comparisons, uncertainty analysis, and biological validation."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "from pathlib import Path\n",
    "import warnings\n",
    "\n",
    "warnings.filterwarnings('ignore')\n",
    "\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import matplotlib.pyplot as plt\n    \n",
    "import seaborn as sns\n",
    "\n",
    "from Bio import SeqIO\n\n",
    "from sklearn.preprocessing import StandardScaler\n",
    "from sklearn.decomposition import PCA\n",
    "from sklearn.cluster import AgglomerativeClustering\n",
    "from sklearn.metrics import silhouette_score\n",
    "\n",
    "print('All libraries loaded successfully.')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "def find_repository_root():\n",
    "    current = Path.cwd().resolve()\n",
    "\n",
    "    candidates = [\n",
    "        current,\n",
    "        current.parent,\n",
    "        current.parent.parent,\n",
    "        current.parent.parent.parent\n",
    "    ]\n\n",
    "    for path in candidates:\n",
    "        if (path / 'data').exists() and (path / 'README.md').exists():\n",
    "            return path\n\n",
    "    raise FileNotFoundError(\n",
    "        'Could not locate the repository root. '\n",
    "        'Run this notebook from the repository or notebooks directory.'\n",
    "    )\n\n",
    "\n",
    "ROOT = find_repository_root()\n\n",
    "DATA_DIR = ROOT / 'data'\n",
    "RAW_DIR = DATA_DIR / 'raw'\n",
    "PROCESSED_DIR = DATA_DIR / 'processed'\n    \n",
    "RESULTS_DIR = ROOT / 'results'\n",
    "FIGURES_DIR = ROOT / 'figures'\n\n",
    "RESULTS_DIR.mkdir(exist_ok=True)\n",
    "FIGURES_DIR.mkdir(exist_ok=True)\n\n",
    "INPUT_FASTA = PROCESSED_DIR / 'TP53_clean.fasta'\n\n",
    "print('Repository:', ROOT)\n",
    "print('Input FASTA:', INPUT_FASTA)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "if not INPUT_FASTA.exists():\n",
    "    raise FileNotFoundError(\n",
    "        f'Input FASTA was not found: {INPUT_FASTA}'\n",
    "    )\n\n",
    "records = list(\n",
    "    SeqIO.parse(INPUT_FASTA, 'fasta')\n",
    ")\n\n",
    "print(f'Sequences loaded: {len(records)}')\n\n",
    "for record in records:\n",
    "    print(f'{record.id}: {len(record.seq)} aa')"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# 🔎 Sequence Quality Control\n",
    "\n",
    "Only standard amino-acid residues are retained for the feature-analysis workflow.\n",
    "\n",
    "This step is intended to make the input data explicit and reproducible before feature extraction."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "VALID_AMINO_ACIDS = set(\n",
    "    'ACDEFGHIKLMNPQRSTVWY'\n",
    ")\n\n",
    "\n",
    "def clean_sequence(sequence):\n",
    "    return (\n",
    "        str(sequence)\n",
    "        .upper()\n",
    "        .replace(' ', '')\n",
    "        .replace('\\n', '')\n",
    "    )\n\n",
    "\n",
    "cleaned_records = []\n\n",
    "for record in records:\n",
    "    sequence = clean_sequence(record.seq)\n    "    \n",
    "    invalid = set(sequence) - VALID_AMINO_ACIDS\n",
    "\n",
    "    if invalid:\n",
    "        print(\n",
    "            f'Excluded {record.id}: '\n",
    "            f'invalid residues = {sorted(invalid)}'\n",
    "        )\n",
    "        continue\n\n",
    "    cleaned_records.append({\n",
    "        'id': record.id,\n",
    "        'description': record.description,\n",
    "        'sequence': sequence\n",
    "    })\n\n",
    "print(\n    'Valid sequences retained:',\n    len(cleaned_records)\n)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "quality_rows = []\n\n",
    "for record in cleaned_records:\n",
    "    sequence = record['sequence']\n",
    "\n",
    "    quality_rows.append({\n",
    "        'id': record['id'],\n    "        'length': len(sequence),\n    "        'unique_amino_acids': len(set(sequence)),\n    "        'unknown_residues': sum(\n    "            aa not in VALID_AMINO_ACIDS\n    "            for aa in sequence\n    "        )\n    "    })\n\n",
    "quality_df = pd.DataFrame(quality_rows)\n\n",
    "quality_df"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# 🧪 Feature Engineering\n",
    "\n",
    "The following sequence-derived features are calculated:\n",
    "\n",
    "- sequence length;\n",
    "- hydrophobic residue fraction;\n",
    "- charged residue fraction;\n",
    "- aromatic residue fraction;\n",
    "- polar residue fraction;\n",
    "- proline fraction;\n",
    "- glycine fraction;\n",
    "- cysteine fraction;\n",
    "- lysine fraction;\n",
    "- arginine fraction;\n",
    "- acidic residue fraction;\n",
    "- basic residue fraction.\n",
    "\n",
    "These features are descriptive sequence properties and should not be interpreted independently as functional or clinical predictors."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "def sequence_features(sequence):\n",
    "    length = len(sequence)\n",
    "\n",
    "    if length == 0:\n",
    "        raise ValueError('Sequence cannot be empty.')\n\n",
    "    hydrophobic = set('AVILMFWY')\n",
    "    charged = set('DEKR')\n",
    "    aromatic = set('FWY')\n",
    "    polar = set('STNQ')\n",
    "    acidic = set('DE')\n",
    "    basic = set('KR')\n\n",
    "    return {\n",
    "        'length': length,\n    \n",
    "        'hydrophobic_fraction': sum(\n",
    "            aa in hydrophobic for aa in sequence\n",
    "        ) / length,\n\n",
    "        'charged_fraction': sum(\n",
    "            aa in charged for aa in sequence\n",
    "        ) / length,\n\n",
    "        'aromatic_fraction': sum(\n",
    "            aa in aromatic for aa in sequence\n",
    "        ) / length,\n\n",
    "        'polar_fraction': sum(\n",
    "            aa in polar for aa in sequence\n",
    "        ) / length,\n\n",
    "        'acidic_fraction': sum(\n",
    "            aa in acidic for aa in sequence\n",
    "        ) / length,\n\n",
    "        'basic_fraction': sum(\n",
    "            aa in basic for aa in sequence\n",
    "        ) / length,\n\n",
    "        'proline_fraction': sequence.count('P') / length,\n    "        'glycine_fraction': sequence.count('G') / length,\n    "        'cysteine_fraction': sequence.count('C') / length,\n    "        'lysine_fraction': sequence.count('K') / length,\n    "        'arginine_fraction': sequence.count('R') / length\n    "    }\n\n",
    "\n",
    "feature_rows = []\n\n",
    "for record in cleaned_records:\n    row = {'id': record['id']}\n    row.update(sequence_features(record['sequence']))\n    feature_rows.append(row)\n\n",
    "feature_df = pd.DataFrame(feature_rows)\n\n",
    "feature_df"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'\n\n",
    "composition_rows = []\n\n",
    "for record in cleaned_records:\n    sequence = record['sequence']\n    length = len(sequence)\n\n",
    "    row = {'id': record['id']}\n\n",
    "    for aa in AMINO_ACIDS:\n",
    "        row[f'aa_{aa}_fraction'] = (\n",
    "            sequence.count(aa) / length\n",
    "        )\n\n",
    "    composition_rows.append(row)\n\n",
    "composition_df = pd.DataFrame(composition_rows)\n\n",
    "composition_df"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "analysis_df = (\n    feature_df\n    .merge(\n        composition_df,\n        on='id',\n        how='left'\n    )\n)\n\n",
    "analysis_df"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# 📊 Exploratory Feature Analysis\n",
    "\n",
    "The next analyses examine relationships among sequence-derived features.\n",
    "\n",
    "These visualizations are exploratory and are not intended to establish causality."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "feature_columns = [\n",
    "    column\n",
    "    for column in feature_df.columns\n",
    "    if column not in ['id']\n",
    "]\n\n",
    "summary_df = (\n",
    "    feature_df[feature_columns]\n",
    "    .describe()\n",
    "    .T\n",
    ")\n\n",
    "summary_df"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "plt.figure(figsize=(12, 7))\n\n",
    "sns.heatmap(\n    feature_df[feature_columns].corr(),\n    annot=True,\n    fmt='.2f',\n    center=0,\n    linewidths=0.5\n)\n\n",
    "plt.title(\n    'Correlation Structure of Sequence-Derived Features'\n)\n\n",
    "plt.tight_layout()\n\n",
    "plt.savefig(\n    FIGURES_DIR / 'eleprotect_feature_correlation.png',\n    dpi=300,\n    bbox_inches='tight'\n)\n\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "plot_features = [\n",
    "    'hydrophobic_fraction',\n",
    "    'charged_fraction',\n",
    "    'aromatic_fraction',\n",
    "    'polar_fraction',\n    "    'acidic_fraction',\n    "    'basic_fraction'\n",
    "]\n\n",
    "plot_data = feature_df.set_index('id')[plot_features]\n\n",
    "plot_data.plot(\n    kind='bar',\n    figsize=(12, 7)\n)\n\n",
    "plt.xlabel('Sequence')\n",
    "plt.ylabel('Fraction')\n",
    "plt.title('Comparative Sequence-Derived Features')\n    \n",
    "plt.xticks(rotation=45, ha='right')\n",
    "plt.tight_layout()\n\n",
    "plt.savefig(\n    FIGURES_DIR / 'eleprotect_sequence_features.png',\n    dpi=300,\n    bbox_inches='tight'\n)\n\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# 🧠 Feature Standardization\n",
    "\n",
    "Features are standardized before dimensionality reduction and clustering.\n",
    "\n",
    "Standardization prevents features with larger numerical scales from dominating the exploratory analysis."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "model_features = feature_df[feature_columns].copy()\n\n",
    "scaler = StandardScaler()\n\n",
    "X_scaled = scaler.fit_transform(model_features)\n\n",
    "print('Scaled feature matrix shape:', X_scaled.shape)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# 📉 Principal Component Analysis\n",
    "\n",
    "Principal Component Analysis (PCA) is used here as an exploratory method to summarize variation across the sequence-derived feature space.\n",
    "\n",
    "The PCA coordinates should be interpreted as mathematical projections rather than biological axes."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "if len(feature_df) >= 2:\n",
    "    n_components = min(\n",
    "        2,\n",
    "        X_scaled.shape[0],\n",
    "        X_scaled.shape[1]\n",
    "    )\n\n",
    "    pca = PCA(\n    "        n_components=n_components,\n    "        random_state=42\n",
    "    )\n\n",
    "    X_pca = pca.fit_transform(X_scaled)\n\n",
    "    print(\n    "        'Explained variance ratio:',\n    "        pca.explained_variance_ratio_\n    "    )\nelse:\n    X_pca = np.empty((len(feature_df), 0))\n    print('PCA requires at least two sequences.')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "if X_pca.shape[1] >= 2:\n    plt.figure(figsize=(9, 7))\n\n",
    "    plt.scatter(\n    "        X_pca[:, 0],\n    "        X_pca[:, 1],\n    "        s=100\n    "    )\n\n",
    "    for i, sequence_id in enumerate(feature_df['id']):\n    "        plt.annotate(\n    "            sequence_id,\n    "            (X_pca[i, 0], X_pca[i, 1]),\n    "            xytext=(5, 5),\n    "            textcoords='offset points',\n    "            fontsize=8\n    "        )\n\n",
    "    plt.xlabel('Principal Component 1')\n    "    plt.ylabel('Principal Component 2')\n    "    plt.title('Exploratory PCA of TP53 Sequence Features')\n\n",
    "    plt.tight_layout()\n\n",
    "    plt.savefig(\n    "        FIGURES_DIR / 'eleprotect_pca.png',\n    "        dpi=300,\n    "        bbox_inches='tight'\n    "    )\n\n",
    "    plt.show()\nelse:\n    print('Two-dimensional PCA visualization is unavailable.')"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# 🔗 Exploratory Clustering\n",
    "\n",
    "Agglomerative clustering is used to examine whether the sequence-derived feature space contains obvious computational groupings.\n",
    "\n",
    "The clustering is exploratory and should not be interpreted as evidence of biological categories or functional classes."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "n_sequences = len(feature_df)\n\n",
    "if n_sequences >= 3:\n    "    n_clusters = min(2, n_sequences - 1)\n\n",
    "    clustering = AgglomerativeClustering(\n    "        n_clusters=n_clusters,\n    "        linkage='ward'\n    "    )\n\n",
    "    cluster_labels = clustering.fit_predict(X_scaled)\n\n",
    "    clustering_df = pd.DataFrame({\n    "        'id': feature_df['id'],\n    "        'cluster': cluster_labels\n    "    })\n\n",
    "    if len(set(cluster_labels)) > 1:\n    "        silhouette = silhouette_score(\n    "            X_scaled,\n    "            cluster_labels\n    "        )\n\n",
    "        print(\n    "            f'Exploratory silhouette score: {silhouette:.4f}'\n    "        )\n    "    else:\n    "        print('Silhouette score unavailable for one cluster.')\nelse:\n    clustering_df = pd.DataFrame({\n    "        'id': feature_df['id'],\n    "        'cluster': np.nan\n    "    })\n\n",
    "    print('Clustering requires at least three sequences.')\n\n",
    "clustering_df"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# 🧮 Exploratory Feature Ranking\n",
    "\n",
    "For visualization purposes, the absolute standardized magnitude of each feature is summarized across the available sequences.\n",
    "\n",
    "This ranking is descriptive only. It does not identify biologically important features or predictive biomarkers."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "feature_importance_proxy = pd.DataFrame({\n    "        'feature': feature_columns,\n    "        'mean_absolute_standardized_value': np.mean(\n    "            np.abs(X_scaled),\n    "            axis=0\n    "        )\n    "    })\n\n",
    "feature_importance_proxy = feature_importance_proxy.sort_values(\n    "    'mean_absolute_standardized_value',\n    "    ascending=False\n    ")\n\n",
    "feature_importance_proxy"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "top_features = feature_importance_proxy.head(10)\n\n",
    "plt.figure(figsize=(10, 7))\n\n",
    "plt.barh(\n    "    top_features['feature'][::-1],\n    "    top_features['mean_absolute_standardized_value'][::-1]\n    ")\n\n",
    "plt.xlabel('Mean absolute standardized magnitude')\n    "plt.ylabel('Feature')\n    "plt.title('Exploratory Sequence Feature Ranking')\n\n",
    "plt.tight_layout()\n\n",
    "plt.savefig(\n    "    FIGURES_DIR / 'eleprotect_feature_ranking.png',\n    "    dpi=300,\n    "    bbox_inches='tight'\n    ")\n\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# 💾 Reproducible Outputs\n",
    "\n",
    "The generated feature table and exploratory clustering results are exported to the repository's `results/` directory.\n\n",
    "Generated figures are stored in `figures/`."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "feature_output = RESULTS_DIR / 'eleprotect_sequence_features.csv'\n\n",
    "feature_df.to_csv(\n    "    feature_output,\n    "    index=False\n    ")\n\n",
    "cluster_output = RESULTS_DIR / 'eleprotect_exploratory_clusters.csv'\n\n",
    "clustering_df.to_csv(\n    "    cluster_output,\n    "    index=False\n    ")\n\n",
    "ranking_output = RESULTS_DIR / 'eleprotect_feature_ranking.csv'\n\n",
    "feature_importance_proxy.to_csv(\n    "    ranking_output,\n    "    index=False\n    ")\n\n",
    "print('Saved outputs:')\n",
    "print(feature_output)\n",
    "print(cluster_output)\n",
    "print(ranking_output)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# ⚠️ Scientific Interpretation and Limitations\n",
    "\n",
    "The analyses in this notebook are exploratory.\n",
    "\n",
    "They should not be interpreted as evidence of:\n",
    "\n",
    "- elephant cancer resistance;\n",
    "- altered TP53 function;\n",
    "- cancer protection;\n",
    "- pathogenicity;\n",
    "- clinical risk;\n",
    "- therapeutic potential;\n",
    "- or biological causality.\n",
    "\n",
    "The available sequence set is small and therefore unsuitable for claiming generalizable machine-learning performance.\n",
    "\n",
    "A future predictive study would require:\n",
    "\n",
    "1. substantially larger datasets;\n",
    "2. independently curated biological labels;\n",
    "3. biologically justified feature selection;\n",
    "4. independent training, validation, and test datasets;\n",
    "5. prevention of data leakage;\n",
    "6. appropriate baseline models;\n",
    "7. external validation;\n",
    "8. uncertainty estimation;\n",
    "9. model interpretability;\n",
    "10. independent biological validation."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# 🚀 Future Research Directions\n",
    "\n",
    "Potential extensions include:\n",
    "\n",
    "- larger cross-species TP53 datasets;\n",
    "- evolutionary conservation features;\n",
    "- structural protein descriptors;\n",
    "- protein language-model embeddings;\n",
    "- supervised functional labels;\n",
    "- interpretable machine-learning models;\n",
    "- external validation datasets;\n",
    "- residue-level attribution;\n",
    "- multi-omic integration;\n",
    "- experimental validation.\n",
    "\n",
    "Any predictive application should remain exploratory until independently validated."
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "name": "python",
   "version": "3.x"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
