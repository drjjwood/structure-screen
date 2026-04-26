from pathlib import Path
import argparse
import subprocess
import pandas as pd


def run(cmd: list[str]) -> None:
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def find_reference_pdb(target_name: str, reference_dir: Path) -> Path | None:
    target_base = target_name.rsplit("_", 1)[0]

    candidates = list(reference_dir.glob("*.pdb"))

    for path in candidates:
        if path.stem == target_name:
            return path

    for path in candidates:
        if path.stem.startswith(target_base) or target_name.startswith(path.stem):
            return path

    return None


def align_structure(
    target_name: str,
    target_pdb: Path,
    query_pdb: Path,
    aligned_dir: Path,
) -> Path:
    aligned_dir.mkdir(parents=True, exist_ok=True)

    out_prefix = aligned_dir / target_name
    aligned_pdb = aligned_dir / f"{target_name}_aligned.pdb"

    run([
        "USalign",
        str(target_pdb),
        str(query_pdb),
        "-o",
        str(out_prefix),
    ])

    produced = Path(f"{out_prefix}.pdb")

    if not produced.exists():
        raise FileNotFoundError(f"Expected aligned file not found: {produced}")

    produced.rename(aligned_pdb)

    return aligned_pdb


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--query-pdb", required=True)
    parser.add_argument("--results-csv", required=True)
    parser.add_argument("--reference-dir", default="inputs/reference_structures")
    parser.add_argument("--aligned-dir", required=True)

    args = parser.parse_args()

    query_pdb = Path(args.query_pdb)
    results_csv = Path(args.results_csv)
    reference_dir = Path(args.reference_dir)
    aligned_dir = Path(args.aligned_dir)

    df = pd.read_csv(results_csv)

    if df.empty:
        print("No hits to align.")
        return

    for target_name in df["target"].unique():
        target_pdb = find_reference_pdb(target_name, reference_dir)

        if target_pdb is None:
            print(f"Could not find PDB for {target_name}")
            continue

        aligned_path = align_structure(
            target_name=target_name,
            target_pdb=target_pdb,
            query_pdb=query_pdb,
            aligned_dir=aligned_dir,
        )

        print(f"Aligned structure written to: {aligned_path}")


if __name__ == "__main__":
    main()