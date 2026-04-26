from pathlib import Path
import subprocess
import shutil


def run_command(cmd: list[str]) -> None:
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def resolve_foldseek_binary() -> str:
    """Find foldseek binary from PATH or known local install path."""
    foldseek_bin = shutil.which("foldseek")
    if foldseek_bin:
        return foldseek_bin

    local_candidate = Path.home() / "tools" / "foldseek" / "foldseek" / "bin" / "foldseek"
    if local_candidate.exists():
        return str(local_candidate)

    raise FileNotFoundError(
        "Could not find 'foldseek'. Add it to PATH or install it at "
        "~/tools/foldseek/foldseek/bin/foldseek."
    )


def main():
    query_pdb = Path("inputs/query.pdb")
    reference_dir = Path("inputs/reference_structures")
    output_dir = Path("outputs")
    tmp_dir = Path("tmp")

    output_dir.mkdir(exist_ok=True)
    tmp_dir.mkdir(exist_ok=True)

    db_path = output_dir / "reference_db"
    results_path = output_dir / "foldseek_results.tsv"

    if not query_pdb.exists():
        raise FileNotFoundError(f"Missing query structure: {query_pdb}")

    if not reference_dir.exists():
        raise FileNotFoundError(f"Missing reference directory: {reference_dir}")

    foldseek = resolve_foldseek_binary()

    # 1. Create Foldseek database from reference structures
    run_command([
        foldseek,
        "createdb",
        str(reference_dir),
        str(db_path),
    ])

    # 2. Search query structure against database
    run_command([
        foldseek,
        "easy-search",
        str(query_pdb),
        str(db_path),
        str(results_path),
        str(tmp_dir),
        "--format-output",
        "query,target,evalue,bits,alntmscore,qtmscore,ttmscore,rmsd,alnlen,qcov,tcov,fident",
    ])

    print(f"\nDone. Results written to: {results_path}")


if __name__ == "__main__":
    main()