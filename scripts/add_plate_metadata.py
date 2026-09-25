#!/usr/bin/env python3
"""Append a plate's rows to data/metadata.csv, deriving cell_id from the FASTQs.

cell_id must equal the FASTQ stem (everything before _R1_001.fastq.gz), which
carries the S-number and lane tag. Those cannot be guessed, so they are read
off the files rather than constructed.

Usage:
  python scripts/add_plate_metadata.py --plate NODPDL1_2 --batch batch5 \
      --strain-group NODPDL1_2 [--dry-run]
"""
import argparse, csv, os, re, sys, collections

# Well -> condition for the NODPDL1_2 sort (see the plate map):
#   A1-A12        CD45+ MHCII+ PDL1+   (12)  <- normalization reference
#   B1-D12, E1-E6 CD45- MHCII+ PDL1+   (42)
#   E7-E12, F-H   CD45- MHCII- PDL1+   (42)
def condition_for(row, col):
    if row == 'A':
        return 'CD45pos_MHCIIpos_PDL1pos'
    if row in 'BCD' or (row == 'E' and col <= 6):
        return 'CD45neg_MHCIIpos_PDL1pos'
    return 'CD45neg_MHCIIneg_PDL1pos'

ap = argparse.ArgumentParser()
ap.add_argument('--plate', required=True)
ap.add_argument('--batch', required=True)
ap.add_argument('--strain-group', required=True)
ap.add_argument('--fastq-dir')
ap.add_argument('--metadata', default='data/metadata.csv')
ap.add_argument('--dry-run', action='store_true')
a = ap.parse_args()

fq = a.fastq_dir
if not fq:
    import yaml
    fq = yaml.safe_load(open('config.yaml'))['FASTQ_DIR']

# Stems for this plate. Matches <PLATE>_<WELL>_..._R1_001.fastq.gz and requires
# the well to follow the plate name exactly, so NODPDL1 does not match
# NODPDL1_2 or vice versa.
pat = re.compile(rf'^{re.escape(a.plate)}_([A-H])(\d{{1,2}})_.*?_R1_001\.fastq\.gz$')
stems = {}
for f in sorted(os.listdir(fq)):
    m = pat.match(f)
    if m:
        well = f"{m.group(1)}{int(m.group(2))}"
        stem = f[:-len('_R1_001.fastq.gz')]
        if well in stems:
            sys.exit(f"ERROR: well {well} matched twice: {stems[well]} and {stem}")
        stems[well] = stem

if not stems:
    sys.exit(f"ERROR: no R1 files matching plate '{a.plate}' in {fq}")

# Every R1 needs its R2, or the pair is unusable downstream.
missing_r2 = [s for s in stems.values()
              if not os.path.exists(os.path.join(fq, s + '_R2_001.fastq.gz'))]
if missing_r2:
    sys.exit(f"ERROR: {len(missing_r2)} R1 file(s) have no R2: {missing_r2[:5]}")

expected = [f"{r}{c}" for r in 'ABCDEFGH' for c in range(1, 13)]
absent = [w for w in expected if w not in stems]
extra  = [w for w in stems if w not in expected]
print(f"{a.plate}: {len(stems)} wells with R1+R2")
if absent: print(f"  WARNING: {len(absent)} well(s) with no FASTQ: {absent}")
if extra:  print(f"  WARNING: unexpected well(s): {extra}")

rows = []
for w in expected:
    if w not in stems: continue
    cond = condition_for(w[0], int(w[1:]))
    rows.append([stems[w], a.plate, w, cond, a.batch, a.strain_group])

print("  condition breakdown:")
for c, n in sorted(collections.Counter(r[3] for r in rows).items()):
    print(f"    {c:32s} {n}")

existing = list(csv.reader(open(a.metadata)))
have = {r[0] for r in existing[1:]}
dupes = [r[0] for r in rows if r[0] in have]
if dupes:
    sys.exit(f"ERROR: {len(dupes)} cell_id(s) already in metadata, e.g. {dupes[:3]}. "
             "Refusing to duplicate rows.")
if any(r[1] == a.plate for r in existing[1:]):
    sys.exit(f"ERROR: plate '{a.plate}' already present in {a.metadata}.")

if a.dry_run:
    print("\n  --dry-run: first 3 and last 3 rows that WOULD be added:")
    for r in rows[:3] + rows[-3:]: print("   ", ",".join(r))
    sys.exit(0)

with open(a.metadata, 'a', newline='') as fh:
    csv.writer(fh).writerows(rows)
print(f"\nAppended {len(rows)} rows to {a.metadata} "
      f"({len(existing)-1} -> {len(existing)-1+len(rows)} cells)")
