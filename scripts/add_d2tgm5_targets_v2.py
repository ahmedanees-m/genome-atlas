
import pandas as pd
from pathlib import Path

DATA = Path.home() / "pen-stack" / "data"
TARGETS_V2 = DATA / "processed/targets_v2.parquet"

df = pd.read_parquet(TARGETS_V2)
print(f"targets_v2.parquet rows before: {len(df)}")
print(f"D2TGM5 present: {'D2TGM5' in df['accession'].values}")

if "D2TGM5" in df["accession"].values:
    print("Already present, skipping")
    exit(0)

sequence = (
    "MEQELHFIGIDVSKAKLDVDVLRPDGRHRSKKFANTPKGHDELLRWLSGHRVAPAHICMEATSTYMEDVAAHLSDAGYTVSVINPALGKAFAQSEGLRSKTD"
    "AVDARMLAEFCRQKRPPAWEAPHPVERALRALVLRHQSLTDMHTQELNRLETAREVQRPSIDAHLLWLHAELKRIEKQIKDLTDDDPDMKHRRKLLESIPGI"
    "GEKTSAVLLAYTGLKERFTHARQFAAFAGLTPRRYESGSSVNRASRMSKAGHASLRRALYMPAMVAVSKTEWGRAFRDRLAGNGKKGKVIIGAMMRKLAQVA"
    "YGVLKSGVPFDASRHNPVAA"
)
new_row = {
    "accession": "D2TGM5",
    "primary_mechanism_bucket": "DSB_FREE_TRANSEST_RECOMBINASE",
    "reviewed": "unreviewed",
    "protein_name": "ISCro4 bridge recombinase (formerly IS622; UniProt D2TGM5)",
    "organism_name": "Citrobacter rodentium ICC168",
    "organism_id": 637910,
    "length": 326,
    "xref_pfam": "PF01548;PF02371;",
    "xref_pdb": None,
    "xref_alphafolddb": "D2TGM5;",
    "sequence": sequence,
    "lineage_ids": "131567 (no rank), 2 (domain), 1783272 (kingdom), 1224 (phylum), 1236 (class), 91347 (order), 543 (family), 544 (genus)",
    "curiosity_score": float("nan"),
}
updated = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
updated.to_parquet(TARGETS_V2, index=False, compression="zstd")
print(f"Added D2TGM5, now {len(updated)} rows")
