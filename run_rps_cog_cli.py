#!/usr/bin/env python3
import os, sys, argparse, subprocess, shutil, csv
from collections import defaultdict

ROOT = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(ROOT))

# Try to locate the repo root (where B-Pan.py is)
for up in [".", "..", "../..", "../../..", "../../../.."]:
    cand = os.path.abspath(os.path.join(ROOT, up, "B-Pan.py"))
    if os.path.exists(cand):
        REPO = os.path.dirname(cand)
        break

TOOLS = [
    os.path.join(REPO, 'Data', 'Tools', 'blast+', 'rpsblast.exe'),
    os.path.join(REPO, 'Data', 'Tools', 'blast-2.17.0+', 'bin', 'rpsblast.exe'),
    'rpsblast'
]

DB_CANDIDATES = [
    os.path.join(REPO, 'CDD_full', 'Cdd'),
    os.path.join(REPO, 'Data', 'Db', 'COG_CDD', 'RPS_COG'),
    os.path.join(REPO, 'Data', 'Db', 'COG_CDD', 'Cdd'),
    os.path.join(REPO, 'CDD', 'Cdd'),
]


def which_rpsblast():
    for t in TOOLS:
        exe = t
        if sys.platform.startswith('win') and not exe.lower().endswith('.exe'):
            exe = exe + '.exe' if not exe.endswith('rpsblast') else exe + '.exe'
        found = shutil.which(exe)
        if found:
            return found
        if os.path.isfile(exe):
            return exe
    return None


def find_rps_db():
    def has_alias(base):
        return os.path.exists(base + '.pal') or any(os.path.exists(base + ext) for ext in ('.pin', '.psq', '.phr', '.rps'))
    for base in DB_CANDIDATES:
        if has_alias(base):
            return base
    return None


def merge_covered_len(intervals):
    if not intervals:
        return 0
    intervals = [(min(a, b), max(a, b)) for a, b in intervals]
    intervals.sort(key=lambda x: x[0])
    merged = []
    s, e = intervals[0]
    for a, b in intervals[1:]:
        if a <= e:
            if b > e:
                e = b
        else:
            merged.append((s, e))
            s, e = a, b
    merged.append((s, e))
    return sum((b - a + 1) for a, b in merged)


def run_rps(query, db_base, out_path, rps_exe):
    outfmt = '6 qseqid sacc pident length qlen qstart qend evalue bitscore'
    cmd = [
        rps_exe,
        '-query', query,
        '-db', db_base,
        '-out', out_path,
        '-outfmt', outfmt
    ]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if res.returncode != 0:
        sys.stderr.write('RPS-BLAST failed (return code %s)\n' % res.returncode)
        sys.stderr.write(res.stderr.decode(errors='ignore'))
    return os.path.exists(out_path) and os.path.getsize(out_path) > 0


def load_cog_type_to_symbol(repo_root):
    path = os.path.join(repo_root, 'Data', 'Db', 'COG_END.csv')
    if not os.path.exists(path):
        return {}
    with open(path, newline='', encoding='utf-8') as f:
        rdr = csv.DictReader(f)
        type_to_symbol = {}
        for row in rdr:
            t = row.get('TYPE')
            s = row.get('SYMBOL')
            if t and s and t not in type_to_symbol:
                type_to_symbol[t] = s
        return type_to_symbol


def parse_and_summarize(tsv_path, id_thr, cov_thr):
    groups = {}
    with open(tsv_path, 'r', encoding='utf-8') as f:
        for ln in f:
            parts = ln.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            q, s = parts[0], parts[1]
            try:
                pid = float(parts[2]); alen = float(parts[3]); qlen = float(parts[4])
                qs = int(parts[5]); qe = int(parts[6])
                evalue = float(parts[7]); bits = float(parts[8])
            except Exception:
                continue
            key = (q, s)
            obj = groups.get(key)
            if obj is None:
                obj = {'qlen': qlen, 'intervals': [], 'pid': pid, 'bits': bits, 'evalue': evalue}
                groups[key] = obj
            else:
                obj['pid'] = max(obj['pid'], pid)
                obj['bits'] = max(obj['bits'], bits)
                obj['evalue'] = min(obj['evalue'], evalue)
                if qlen > 0:
                    obj['qlen'] = qlen
            obj['intervals'].append((qs, qe))

    passed = []
    for (q, s), g in groups.items():
        covered = merge_covered_len(g['intervals'])
        qlen = g['qlen'] if g['qlen'] else 0.0
        cov = (covered / qlen) * 100.0 if qlen > 0 else 0.0
        if g['pid'] >= id_thr and cov >= cov_thr:
            passed.append({'QUERY': q, 'SUBJECT': s, 'PERCENTAGE IDENTITY': g['pid'], 'BIT SCORE': g['bits'], 'E VALUE': g['evalue'], 'COVERAGE': cov})
    return passed


def write_outputs(repo_root, rows, type_to_symbol):
    out_dir = os.path.join(repo_root, 'Results', 'COG')
    os.makedirs(out_dir, exist_ok=True)
    sym_csv = os.path.join(out_dir, 'COG_SYMBOL.csv')
    perc_csv = os.path.join(out_dir, 'COG_PERC.csv')

    # Map symbols
    for r in rows:
        r['SYMBOL'] = type_to_symbol.get(r['SUBJECT'], 'Unknown')

    # Write symbol assignments
    fields = ['QUERY', 'SUBJECT', 'SYMBOL', 'PERCENTAGE IDENTITY', 'BIT SCORE', 'E VALUE', 'COVERAGE']
    with open(sym_csv, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)

    # Summary counts
    counts = defaultdict(int)
    for r in rows:
        counts[r['SYMBOL']] += 1
    total = sum(counts.values())
    with open(perc_csv, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['Functional Category', 'Number of Sequences', 'Percentage'])
        for sym, n in sorted(counts.items(), key=lambda x: -x[1]):
            pct = (n / total * 100.0) if total else 0.0
            w.writerow([sym, n, f"{pct:.2f}"])

    return sym_csv, perc_csv


def main():
    ap = argparse.ArgumentParser(description='Run RPS-BLAST COG annotation (CLI)')
    ap.add_argument('--query', required=True)
    ap.add_argument('--db', default=None, help='RPS DB base; if omitted, autodetect')
    ap.add_argument('--id', type=float, default=40.0)
    ap.add_argument('--cov', type=float, default=50.0)
    args = ap.parse_args()

    repo = REPO
    rps_exe = which_rpsblast()
    if not rps_exe:
        print('ERROR: rpsblast.exe not found. Expected under Data/Tools/blast-2.17.0+/bin or on PATH.', file=sys.stderr)
        return 2
    db_base = args.db or find_rps_db()
    if not db_base:
        print('ERROR: RPS database not found. Run Data/Db/COG_CDD/setup_cog_rps_db.py first.', file=sys.stderr)
        return 3

    out_dir = os.path.join(repo, 'Results', 'COG')
    os.makedirs(out_dir, exist_ok=True)
    tsv = os.path.join(out_dir, 'rps_cli_out.tsv')

    ok = run_rps(args.query, db_base, tsv, rps_exe)
    if not ok:
        print('ERROR: RPS-BLAST produced no output.', file=sys.stderr)
        return 4

    rows = parse_and_summarize(tsv, args.id, args.cov)
    if not rows:
        print('No hits passed thresholds (identity>=%.1f, coverage>=%.1f).' % (args.id, args.cov))
        return 0

    type_to_symbol = load_cog_type_to_symbol(repo)
    sym_csv, perc_csv = write_outputs(repo, rows, type_to_symbol)
    print('Wrote:', sym_csv)
    print('Wrote:', perc_csv)
    return 0

if __name__ == '__main__':
    sys.exit(main())

