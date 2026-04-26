from pathlib import Path
import argparse
import shutil
import subprocess
import sys


REFERENCE_DIR = Path("inputs/reference_structures")
OUTPUTS_DIR = Path("outputs")
TMP_DIR = Path("tmp")


def run(cmd: list[str]) -> None:
    print("\nRunning:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def resolve_binary(name: str, fallback: Path | None = None) -> str:
    """Resolve a binary path from PATH, with optional fallback."""
    bin_path = shutil.which(name)
    if bin_path:
        return bin_path

    if fallback is not None and fallback.exists():
        return str(fallback)

    if fallback is not None:
        raise FileNotFoundError(
            f"Could not find '{name}' on PATH or at: {fallback}"
        )

    raise FileNotFoundError(f"Could not find '{name}' on PATH")


def parse_fasta(fasta_path: Path) -> dict[str, str]:
    records = {}
    current_id = None
    current_seq = []

    with open(fasta_path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(current_seq)

                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        records[current_id] = "".join(current_seq)

    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}")

    return records


def write_single_fasta(query_id: str, sequence: str, out_path: Path) -> None:
    with open(out_path, "w") as f:
        f.write(f">{query_id}\n")
        f.write(f"{sequence}\n")


def build_reference_db(reference_db: Path) -> None:
    if not REFERENCE_DIR.exists():
        raise FileNotFoundError(f"Missing reference structure directory: {REFERENCE_DIR}")

    if not any(REFERENCE_DIR.glob("*.pdb")):
        raise FileNotFoundError(f"No .pdb files found in {REFERENCE_DIR}")

    foldseek = resolve_binary(
        "foldseek",
        Path.home() / "tools" / "foldseek" / "foldseek" / "bin" / "foldseek",
    )

    run([
        foldseek,
        "createdb",
        str(REFERENCE_DIR),
        str(reference_db),
    ])


def run_colabfold(query_fasta: Path, query_dir: Path) -> Path:
    colabfold_dir = query_dir / "colabfold"
    colabfold_dir.mkdir(parents=True, exist_ok=True)

    # Prefer LocalColabFold's binary (known-good complete install) over any
    # potentially incomplete colabfold_batch on PATH.
    local_colabfold_batch = (
        Path.home() / "localcolabfold" / ".pixi" / "envs" / "default" / "bin" / "colabfold_batch"
    )
    if local_colabfold_batch.exists():
        colabfold_batch = str(local_colabfold_batch)
    else:
        colabfold_batch = resolve_binary("colabfold_batch")

    run([
        colabfold_batch,
        str(query_fasta),
        str(colabfold_dir),
    ])

    pdb_candidates = sorted(colabfold_dir.glob("*relaxed_rank_001*.pdb"))

    if not pdb_candidates:
        pdb_candidates = sorted(colabfold_dir.glob("*rank_001*.pdb"))

    if not pdb_candidates:
        raise FileNotFoundError(f"No ColabFold PDB output found in {colabfold_dir}")

    query_pdb = query_dir / "query.pdb"
    shutil.copyfile(pdb_candidates[0], query_pdb)

    return query_pdb


def copy_existing_query_pdb(query_id: str, query_dir: Path, existing_pdb_dir: Path) -> Path:
    """
    For hackathon/testing mode.

    Looks for:
      existing_pdb_dir/<query_id>.pdb

    Example:
      inputs/query_pdbs/protein1.pdb
    """
    source_pdb = existing_pdb_dir / f"{query_id}.pdb"

    if not source_pdb.exists():
        raise FileNotFoundError(
            f"Expected existing PDB for {query_id}: {source_pdb}\n"
            "Either add this file or run with --run-colabfold."
        )

    query_pdb = query_dir / "query.pdb"
    shutil.copyfile(source_pdb, query_pdb)

    return query_pdb


def run_foldseek_search(query_pdb: Path, reference_db: Path, query_dir: Path) -> Path:
    results_path = query_dir / "foldseek_results.tsv"
    query_tmp = TMP_DIR / query_dir.name
    query_tmp.mkdir(parents=True, exist_ok=True)

    foldseek = resolve_binary(
        "foldseek",
        Path.home() / "tools" / "foldseek" / "foldseek" / "bin" / "foldseek",
    )

    run([
        foldseek,
        "easy-search",
        str(query_pdb),
        str(reference_db),
        str(results_path),
        str(query_tmp),
        "--format-output",
        "query,target,evalue,bits,alntmscore,qtmscore,ttmscore,rmsd,alnlen,qcov,tcov,fident",
    ])

    return results_path


def score_hits(query_dir: Path) -> Path:
    run([
        sys.executable,
        "scripts/score_foldseek_hits.py",
        "--results",
        str(query_dir / "foldseek_results.tsv"),
        "--out-csv",
        str(query_dir / "scored_similarity_hits.csv"),
        "--out-report",
        str(query_dir / "similarity_report.txt"),
    ])

    return query_dir / "scored_similarity_hits.csv"


def align_hits(query_dir: Path) -> None:
    scored_csv = query_dir / "scored_similarity_hits.csv"

    if not scored_csv.exists():
        print(f"No scored hits found for {query_dir.name}; skipping alignment.")
        return

    run([
        sys.executable,
        "scripts/align_hit_structures.py",
        "--query-pdb",
        str(query_dir / "query.pdb"),
        "--results-csv",
        str(scored_csv),
        "--reference-dir",
        str(REFERENCE_DIR),
        "--aligned-dir",
        str(query_dir / "aligned_structures"),
    ])


def launch_dashboard() -> None:
    run([
        "streamlit",
        "run",
        "scripts/app.py",
    ])


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta",
        required=True,
        help="Input FASTA containing one or more protein sequences.",
    )

    parser.add_argument(
        "--run-colabfold",
        action="store_true",
        help="Run ColabFold to generate query structures.",
    )

    parser.add_argument(
        "--existing-pdb-dir",
        default="inputs/query_pdbs",
        help="Directory containing precomputed query PDBs named <query_id>.pdb.",
    )

    parser.add_argument(
        "--skip-dashboard",
        action="store_true",
        help="Run pipeline but do not launch Streamlit.",
    )

    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    existing_pdb_dir = Path(args.existing_pdb_dir)

    OUTPUTS_DIR.mkdir(exist_ok=True)
    TMP_DIR.mkdir(exist_ok=True)

    reference_db = OUTPUTS_DIR / "reference_db"

    print("Building Foldseek reference database...")
    build_reference_db(reference_db)

    records = parse_fasta(fasta_path)

    for query_id, sequence in records.items():
        print(f"\n==============================")
        print(f"Processing query: {query_id}")
        print(f"==============================")

        query_dir = OUTPUTS_DIR / query_id
        query_dir.mkdir(parents=True, exist_ok=True)

        query_fasta = query_dir / f"{query_id}.fasta"
        write_single_fasta(query_id, sequence, query_fasta)

        if args.run_colabfold:
            query_pdb = run_colabfold(query_fasta, query_dir)
        else:
            query_pdb = copy_existing_query_pdb(
                query_id=query_id,
                query_dir=query_dir,
                existing_pdb_dir=existing_pdb_dir,
            )

        run_foldseek_search(query_pdb, reference_db, query_dir)
        score_hits(query_dir)
        align_hits(query_dir)

    print("\nPipeline complete.")

    if not args.skip_dashboard:
        launch_dashboard()


if __name__ == "__main__":
    main()