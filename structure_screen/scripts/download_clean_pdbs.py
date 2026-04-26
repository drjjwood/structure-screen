from pathlib import Path
import requests

PDBS = {
    "query_pdbs": {
        "1QO2": "query2",  # HisA / ProFAR isomerase
    },
    "reference_structures": {
        "1THF": "hisf_tim_barrel_sequence_diverse_control",
        "2CDS": "lysozyme_negative_control",
    },
}

# PDBS = {
#     "query_pdbs": {
#         # Query for the structure-diverse control
#         "4Y8F": "protein1",
#     },
#     "reference_structures": {
#         # TIM-barrel / related fold controls
#         "1J5T": "indole_3_glycerol_phosphate_synthase_tim_barrel",
#         "1PII": "pra_isomerase_beta_alpha_barrel",

#         # Negative control
#         "2CDS": "lysozyme_negative_control",

#         # Optional old serine protease controls
#         "2CHA": "chymotrypsin_close_homolog",
#         "1EST": "elastase_related_protease",
#     },
# }

# PDBS = {
#     "query": {
#         "1PTN": "trypsin",
#     },
#     "reference_structures": {
#         "2CHA": "chymotrypsin_close_homolog",
#         "1EST": "elastase_related_protease",
#         "2CDS": "lysozyme_negative_control",
#     },
# }


def download_pdb(pdb_id: str, out_path: Path) -> None:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url, timeout=30)

    if response.status_code != 200:
        raise RuntimeError(f"Failed to download {pdb_id} from {url}")

    out_path.write_text(response.text)


def keep_protein_atoms_only(
    in_path: Path,
    out_path: Path,
    chain_id: str | None = None,
) -> None:
    cleaned_lines = []

    with open(in_path) as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue

            if chain_id is not None:
                pdb_chain = line[21].strip()
                if pdb_chain != chain_id:
                    continue

            cleaned_lines.append(line)

    cleaned_lines.append("END\n")
    out_path.write_text("".join(cleaned_lines))


def main():
    base_dir = Path("inputs")
    raw_dir = base_dir / "raw_pdbs"
    query_pdb_dir = base_dir / "query_pdbs"
    ref_dir = base_dir / "reference_structures"

    raw_dir.mkdir(parents=True, exist_ok=True)
    query_pdb_dir.mkdir(parents=True, exist_ok=True)
    ref_dir.mkdir(parents=True, exist_ok=True)

    for group, pdbs in PDBS.items():
        for pdb_id, name in pdbs.items():
            raw_path = raw_dir / f"{pdb_id}_{name}.pdb"

            print(f"Downloading {pdb_id} ({name})...")
            download_pdb(pdb_id, raw_path)

            if group == "query_pdbs":
                clean_path = query_pdb_dir / f"{name}.pdb"
            else:
                clean_path = ref_dir / f"{pdb_id}_{name}.pdb"

            print(f"Cleaning {raw_path} -> {clean_path}")
            keep_protein_atoms_only(raw_path, clean_path)

    print("\nDone.")
    print("Created:")
    print("  inputs/query_pdbs/protein1.pdb")
    print("  inputs/reference_structures/*.pdb")


if __name__ == "__main__":
    main()