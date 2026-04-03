# app.py

import streamlit as st
import pandas as pd
import io
import requests
from Bio import Entrez, SeqIO
import plotly.graph_objects as go

from utils import (
    parse_upload,
    parse_fasta_text,
    clean_sequence,
    translate_if_needed,
    align_and_map,
    build_model_features,
    annotate_hotspots
)

from model import predict_score

# =====================================================
# CONFIG
# =====================================================

st.set_page_config(page_title="EleProtect", page_icon="🧬", layout="wide")

st.markdown("""
<style>
.main { background-color: #0f172a; }
section[data-testid="stSidebar"] { background-color: #111827; }
h1, h2, h3 { color: #38bdf8; }
.stButton>button {
    background: linear-gradient(90deg,#0ea5e9,#22d3ee);
    color: black;
    border-radius: 10px;
    font-weight: bold;
}
.stDownloadButton>button {
    background-color: #10b981;
    color: white;
    border-radius: 10px;
}
</style>
""", unsafe_allow_html=True)

Entrez.email = "your_email@example.com"

# =====================================================
# DEFAULT HUMAN TP53
# =====================================================

def DEFAULT_HUMAN_SEQ():
    return (
        "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
        "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTA"
        "KSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHH"
        "ERCSDPSDGSLAPPQHLIRVEGNLRAEYLDDSITLRHSVVVPYEPPEVGSDCTTIHYNYM"
        "CNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHH"
        "ELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKE"
        "PGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
    )

# =====================================================
# HEADER
# =====================================================

col_logo, col_title = st.columns([1, 6])

with col_logo:
    st.image("EleProtect.png", width=100)

with col_title:
    st.title("EleProtect")
    st.caption("Comparative TP53 Hotspot & Evolutionary Analysis Platform")

st.markdown("---")

# =====================================================
# DATABASE FETCH
# =====================================================

st.subheader("🔎 Fetch from Database")

col1, col2, col3 = st.columns([2, 3, 1])

with col1:
    db_type = st.selectbox("Database", ["UniProt", "GenBank (NCBI)"])

with col2:
    db_id = st.text_input("Accession ID (e.g., P04637)")

with col3:
    fetch_btn = st.button("Fetch")

fetched_sequence = None

if fetch_btn and db_id:
    try:
        if db_type == "GenBank (NCBI)":
            handle = Entrez.efetch(db="protein", id=db_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            fetched_sequence = str(record.seq)
            st.success("Sequence fetched from GenBank")

        elif db_type == "UniProt":
            fasta_url = f"https://rest.uniprot.org/uniprotkb/{db_id}.fasta"
            r = requests.get(fasta_url)
            if r.status_code == 200:
                record = list(SeqIO.parse(io.StringIO(r.text), "fasta"))[0]
                fetched_sequence = str(record.seq)
                st.success("Sequence fetched from UniProt")

        if fetched_sequence:
            st.metric("Protein Length", f"{len(fetched_sequence)} aa")

    except Exception as e:
        st.error(f"Fetch error: {e}")

# =====================================================
# INPUT
# =====================================================

st.markdown("---")
st.subheader("📂 Upload or Paste Sequence")

uploaded_file = st.file_uploader(
    "Upload FASTA / TXT / DOCX / XLSX",
    type=["fasta", "fa", "txt", "docx", "xls", "xlsx"]
)

pasted = st.text_area("Or paste sequence")

hot_input = st.text_input(
    "Hotspot Codons (comma separated)",
    value="175,245,248,249,273,282"
)

run_btn = st.button("Run Analysis")

try:
    hotspots = [int(x.strip()) for x in hot_input.split(",")]
except:
    hotspots = [175,245,248,249,273,282]

# =====================================================
# DOMAIN VISUALIZATION FUNCTION
# =====================================================

def plot_tp53_domains(seq_length):

    fig = go.Figure()

    domains = [
        ("Transactivation", 1, 60),
        ("DNA-binding", 100, 300),
        ("Tetramerization", 320, 360)
    ]

    for name, start, end in domains:
        fig.add_trace(go.Bar(
            x=[end-start],
            y=[name],
            orientation="h",
            base=start
        ))

    fig.update_layout(
        title="TP53 Domain Architecture",
        xaxis_title="Amino Acid Position",
        showlegend=False
    )

    return fig

# =====================================================
# MAIN ANALYSIS
# =====================================================

if run_btn:

    sequences = {}

    if fetched_sequence:
        sequences = {"Fetched": fetched_sequence}
    elif uploaded_file:
        sequences = parse_upload(uploaded_file)
    elif pasted:
        if ">" in pasted:
            sequences = parse_fasta_text(pasted)
        else:
            sequences = {"Manual": pasted}
    else:
        st.error("Provide a sequence first.")
        st.stop()

    results = []

    for name, seq in sequences.items():

        seq_clean = clean_sequence(seq)
        seq_clean = translate_if_needed(seq_clean)

        df_map = align_and_map(
            human_seq=DEFAULT_HUMAN_SEQ(),
            query_seq=seq_clean,
            hotspots=hotspots
        )

        df_map = annotate_hotspots(df_map)

        # Add mutation highlight
        if not df_map.empty:
            df_map["Status"] = df_map.apply(
                lambda x: "⚠ Mutation" if not x["Conserved"] else "✔ Conserved",
                axis=1
            )

        features = build_model_features(df_map)
        score = predict_score(pd.DataFrame([features]))

        results.append({
            "Sequence": name,
            "Conservation_%": round(features["conservation_score"]*100,2),
            "ML_Score": round(float(score),4),
            "Mapping": df_map,
            "Length": len(seq_clean)
        })

    # =====================================================
    # RESULTS TABS
    # =====================================================

    st.markdown("---")

    tabs = st.tabs(["📊 Summary", "🧬 Hotspot Mapping", "📈 Domain Architecture"])

    # SUMMARY TAB
    with tabs[0]:
        summary_df = pd.DataFrame([
            {
                "Sequence": r["Sequence"],
                "Conservation_%": r["Conservation_%"],
                "ML_Score": r["ML_Score"]
            }
            for r in results
        ])
        st.dataframe(summary_df, use_container_width=True)

        avg_cons = summary_df["Conservation_%"].mean()

        if avg_cons > 80:
            st.success("High evolutionary conservation detected.")
        elif avg_cons > 50:
            st.warning("Moderate conservation observed.")
        else:
            st.error("Low conservation — functional alteration possible.")

    # HOTSPOT TAB
    with tabs[1]:
        for r in results:
            st.subheader(r["Sequence"])
            st.dataframe(r["Mapping"], use_container_width=True)

    # DOMAIN TAB
    with tabs[2]:
        for r in results:
            st.subheader(r["Sequence"])
            st.plotly_chart(
                plot_tp53_domains(r["Length"]),
                use_container_width=True
            )

    # Excel Export
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="Summary", index=False)
        for r in results:
            r["Mapping"].to_excel(writer, sheet_name=r["Sequence"][:25], index=False)

    buffer.seek(0)

    st.download_button(
        "Download Excel Report",
        buffer.read(),
        "eleprotect_results.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )

st.markdown("---")
st.caption("EleProtect | Comparative TP53 Genomics Platform")