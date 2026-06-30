# BPAN

BPAN is a pan-genome analysis tool with functional annotation (COG) for bacterial proteomes.

> **Note:** Currently developed and tested for **Windows only**.

## Setup

1. Download `BPAN_Data.zip` from: **https://zenodo.org/records/21079899**
2. Extract it so the contents sit directly inside the repo as `Data/`:
   ```
   BPAN/
   ├── bpan.py
   ├── Data/          <-- extracted here
   │   ├── Db/
   │   ├── Tools/
   │   └── Images/
   ├── Results/
   └── TestData/
   ```


## Project structure

| Path | Description |
|---|---|
| `bpan.py` | Main application |
| `run_rps_cog_cli.py` | CLI helper for RPS-BLAST/COG |
| `Data/Db/` | Reference databases (COG, CDD) + small annotation CSVs |
| `Data/Tools/` | Bundled Windows binaries (DIAMOND, BLAST+, Clustal Omega, MAFFT) |
| `Data/Images/` | UI assets |
| `Results/` | Analysis output (created at runtime) |
| `TestData/` | Sample input files |

## Features

- Pan-genome analysis (core, accessory, unique genes)
- COG functional annotation (RPS-BLAST profile-based and DIAMOND-based)
- Phylogenetic analysis (MSA via Clustal Omega / MAFFT, tree building)


