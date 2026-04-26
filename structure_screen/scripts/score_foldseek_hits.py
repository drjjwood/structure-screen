from pathlib import Path
import argparse
import pandas as pd


COLUMNS = [
    "query",
    "target",
    "evalue",
    "bits",
    "alntmscore",
    "qtmscore",
    "ttmscore",
    "rmsd",
    "alnlen",
    "qcov",
    "tcov",
    "fident",
]


def classify_hit(row) -> str:
    if structural_remote_hit(row):
        return "structure_similar_sequence_diverse"

    if row["qtmscore"] >= 0.7 and row["qcov"] >= 0.7 and row["evalue"] <= 1e-5:
        return "high_similarity"

    if row["qtmscore"] >= 0.5 and row["qcov"] >= 0.5 and row["evalue"] <= 1e-3:
        return "moderate_similarity"

    if row["qtmscore"] >= 0.3:
        return "weak_similarity"

    return "low_similarity"

def structural_remote_hit(row) -> bool:
    return (
        row["qtmscore"] >= 0.5
        and row["qcov"] >= 0.5
        and row["fident"] <= 0.30
        and row["evalue"] <= 1e-3
    )

def explain_hit(row) -> str:
    return (
        f"Target {row['target']} has query-normalised TM-score "
        f"{row['qtmscore']:.3f}, query coverage {row['qcov']:.3f}, "
        f"target coverage {row['tcov']:.3f}, sequence identity "
        f"{row['fident']:.3f}, RMSD {row['rmsd']:.2f}, and E-value "
        f"{row['evalue']:.2e}. Classified as {row['similarity_class']}."
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--results", default="outputs/foldseek_results.tsv")
    parser.add_argument("--out-csv", default="outputs/scored_similarity_hits.csv")
    parser.add_argument("--out-report", default="outputs/similarity_report.txt")

    args = parser.parse_args()

    results_path = Path(args.results)
    scored_path = Path(args.out_csv)
    report_path = Path(args.out_report)

    if not results_path.exists():
        raise FileNotFoundError(f"No Foldseek results found: {results_path}")

    df = pd.read_csv(results_path, sep="\t", names=COLUMNS)

    if df.empty:
        scored_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(scored_path, index=False)
        report_path.write_text("No Foldseek hits found.\n")
        return

    numeric_cols = [
        "evalue",
        "bits",
        "alntmscore",
        "qtmscore",
        "ttmscore",
        "rmsd",
        "alnlen",
        "qcov",
        "tcov",
        "fident",
    ]

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["similarity_class"] = df.apply(classify_hit, axis=1)
    df["explanation"] = df.apply(explain_hit, axis=1)

    df = df.sort_values(
        by=["qtmscore", "qcov", "bits"],
        ascending=False,
    )

    scored_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(scored_path, index=False)

    with open(report_path, "w") as f:
        f.write("Structural similarity report\n")
        f.write("============================\n\n")

        for _, row in df.head(20).iterrows():
            f.write(row["explanation"] + "\n\n")

    print(f"Scored table written to: {scored_path}")
    print(f"Text report written to: {report_path}")


if __name__ == "__main__":
    main()