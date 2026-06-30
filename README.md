# BPAN

BPAN is a pan-genome analysis tool with functional annotation (COG) for bacterial proteomes.

> **Note:** Currently developed and tested for **Windows only**.

## Setup

### 1. Clone the repository
```powershell
git clone https://github.com/<your-username>/BPAN.git
cd BPAN
```

### 2. Install Python dependencies
```powershell
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

### 3. Download the Data folder (databases + bundled tools)
The `Data/` folder (~3+ GB: COG databases, RPS-BLAST profiles, bundled DIAMOND/BLAST+/Clustal
Omega/MAFFT binaries) is hosted separately on Zenodo since it's too large for GitHub.

1. Download `BPAN_Data.zip` from: **[Zenodo link goes here]**
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

### 4. Run BPAN
```powershell
python bpan.py
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

## License

[Add your chosen license here]
