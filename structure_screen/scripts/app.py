from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
import plotly.express as px
import py3Dmol
from stmol import showmol

# RESULTS_PATH = Path("outputs/scored_similarity_hits.csv")
# QUERY_PDB = Path("inputs/query.pdb")
# REFERENCE_DIR = Path("inputs/reference_structures")
# ALIGNED_DIR = Path("outputs/aligned_structures")

BASE_OUTPUTS_DIR = Path("outputs")
REFERENCE_DIR = Path("inputs/reference_structures")

def get_query_dirs() -> list[Path]:
    return sorted(
        [
            path for path in BASE_OUTPUTS_DIR.iterdir()
            if path.is_dir() and (path / "scored_similarity_hits.csv").exists()
        ]
    )


query_dirs = get_query_dirs()

if not query_dirs:
    st.error("No query result folders found in outputs/. Run run_pipeline.py first.")
    st.stop()

selected_query_dir = st.sidebar.selectbox(
    "Choose query",
    query_dirs,
    format_func=lambda p: p.name,
)

RESULTS_PATH = selected_query_dir / "scored_similarity_hits.csv"
QUERY_PDB = selected_query_dir / "query.pdb"
ALIGNED_DIR = selected_query_dir / "aligned_structures"


st.set_page_config(
    page_title="Structure similarity dashboard",
    layout="wide",
)


def load_results():
    df = pd.read_csv(RESULTS_PATH)

    priority_order = {
        "structure_similar_sequence_diverse": 4,
        "high_similarity": 3,
        "moderate_similarity": 2,
        "weak_similarity": 1,
        "low_similarity": 0,
    }

    df["priority"] = df["similarity_class"].map(priority_order).fillna(0)
    df = df.sort_values(
        ["priority", "qtmscore", "qcov"],
        ascending=False,
    )

    return df

def get_binary_decision(df: pd.DataFrame) -> tuple[str, str]:
    flag_hits = df[
        (df["qtmscore"] >= 0.5)
        & (df["qcov"] >= 0.5)
        & (df["fident"] <= 0.30)
        & (df["evalue"] <= 1e-3)
    ]

    if len(flag_hits) > 0:
        return (
            "FLAG",
            f"{len(flag_hits)} hit(s) meet: TM-score ≥ 0.5, qcov ≥ 0.5, sequence identity ≤ 0.30.",
        )

    return (
        "PASS",
        "No hits meet the structure-similar / sequence-diverse threshold.",
    )


def show_binary_decision(decision: str, reason: str) -> None:
    if decision == "FLAG":
        st.error(f"🚩 **FLAG** — {reason}")
    else:
        st.success(f"✅ **PASS** — {reason}")

def aligned_overlay_viewer(query_path: Path, aligned_target_path: Path):
    st.subheader("Aligned structural overlay")

    if aligned_target_path is None or not aligned_target_path.exists():
        st.warning("Aligned structure not found. Run scripts/align_hit_structures.py first.")
        return

    query_text = query_path.read_text()
    target_text = aligned_target_path.read_text()

    view = py3Dmol.view(width=900, height=600)

    view.addModel(query_text, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "blue"}})

    view.addModel(target_text, "pdb")
    view.setStyle({"model": 1}, {"cartoon": {"color": "orange"}})

    view.zoomTo()
    components.html(view._make_html(), height=600, width=900)

def find_aligned_pdb(target_name: str) -> Path | None:
    path = ALIGNED_DIR / f"{target_name}_aligned.pdb"
    return path if path.exists() else None

def find_reference_pdb(target_name: str) -> Path | None:
    """
    Tries to match Foldseek target name back to a PDB file.
    Handles targets like:
      2CHA_chymotrypsin_close_homolog_B
    by matching the start of the filename.
    """
    target_base = target_name.rsplit("_", 1)[0]

    candidates = list(REFERENCE_DIR.glob("*.pdb"))

    for path in candidates:
        if path.stem == target_name:
            return path

    for path in candidates:
        if path.stem.startswith(target_base) or target_name.startswith(path.stem):
            return path

    return None


def structure_viewer(pdb_path: Path, title: str):
    st.subheader(title)

    if pdb_path is None or not pdb_path.exists():
        st.warning("PDB file not found.")
        return

    pdb_text = pdb_path.read_text()

    view = py3Dmol.view(width=700, height=500)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.zoomTo()

    showmol(view, height=500, width=700)


def two_structure_viewer(query_path: Path, target_path: Path):
    """
    Shows query and reference in the same viewer.

    Important: this does not structurally superpose them by itself.
    For a true overlay, you need aligned/superposed PDB coordinates
    from TM-align, PyMOL, US-align, or another alignment tool.
    """
    st.subheader("Query and reference in one viewer")

    if not query_path.exists() or target_path is None or not target_path.exists():
        st.warning("Could not find both PDB files.")
        return

    query_text = query_path.read_text()
    target_text = target_path.read_text()

    view = py3Dmol.view(width=900, height=600)
    view.addModel(query_text, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "blue"}})

    view.addModel(target_text, "pdb")
    view.setStyle({"model": 1}, {"cartoon": {"color": "orange"}})

    view.zoomTo()
    showmol(view, height=600, width=900)


df = load_results()

st.title("Structure Similarity Dashboard")
st.caption(f"Showing results for: `{selected_query_dir.name}`")

decision, reason = get_binary_decision(df)

st.markdown("## Structural Screening decision")
show_binary_decision(decision, reason)

structure_diverse = (
    df["similarity_class"] == "structure_similar_sequence_diverse"
).sum()
high = (df["similarity_class"] == "high_similarity").sum()
moderate = (df["similarity_class"] == "moderate_similarity").sum()
weak = (df["similarity_class"] == "weak_similarity").sum()
low = (df["similarity_class"] == "low_similarity").sum()

col1, col2, col3, col4, col5 = st.columns(5)
col1.metric("Total reported hits", len(df))
col2.metric("Structure-similar / sequence-diverse", structure_diverse)
col3.metric("High similarity", high)
col4.metric("Moderate similarity", moderate)
col5.metric("Low / weak similarity", low + weak)

st.markdown("## Hit overview")

fig = px.scatter(
    df,
    x="fident",
    y="qtmscore",
    size="qcov",
    color="similarity_class",
    hover_name="target",
    hover_data=["evalue", "rmsd", "tcov", "alnlen"],
    labels={
        "fident": "Sequence identity",
        "qtmscore": "Query-normalised TM-score",
        "qcov": "Query coverage",
    },
)

st.plotly_chart(fig, use_container_width=True)

st.dataframe(
    df[
        [
            "target",
            "similarity_class",
            "qtmscore",
            "qcov",
            "tcov",
            "fident",
            "rmsd",
            "evalue",
        ]
    ],
    use_container_width=True,
)

st.markdown("## Inspect individual hit")

target = st.selectbox("Choose a hit", df["target"].tolist())
row = df[df["target"] == target].iloc[0]

st.markdown("### Interpretation")

st.write(row["explanation"])

metric_cols = st.columns(6)
metric_cols[0].metric("TM-score", f"{row['qtmscore']:.3f}")
metric_cols[1].metric("Query coverage", f"{row['qcov']:.3f}")
metric_cols[2].metric("Target coverage", f"{row['tcov']:.3f}")
metric_cols[3].metric("Seq identity", f"{row['fident']:.3f}")
metric_cols[4].metric("RMSD", f"{row['rmsd']:.2f} Å")
metric_cols[5].metric("E-value", f"{row['evalue']:.2e}")

with st.expander("What do these metrics mean?"):
    st.markdown(
        """
        **TM-score** measures structural similarity. Values above ~0.5 often indicate
        a shared fold; values above ~0.7 are strong structural matches.

        **Query coverage** is the fraction of the query structure included in the alignment.

        **Target coverage** is the fraction of the reference structure included in the alignment.

        **Sequence identity** is the fraction of identical residues in the aligned region.

        **RMSD** measures the average atomic distance after alignment. Lower is better,
        but it should not be interpreted without coverage/alignment length.

        **E-value** estimates how likely a match of this strength is by chance. Lower is better.
        """
    )

target_pdb = find_reference_pdb(target)

left, right = st.columns(2)
with left:
    structure_viewer(QUERY_PDB, "Query structure")

with right:
    structure_viewer(target_pdb, "Reference structure")

aligned_pdb = find_aligned_pdb(target)
aligned_overlay_viewer(QUERY_PDB, aligned_pdb)