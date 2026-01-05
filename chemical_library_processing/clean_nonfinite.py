#!/usr/bin/env python3
import argparse, math, csv, sys
from rdkit import Chem
from rdkit.Chem import AllChem

def has_finite_coords(mol) -> bool:
    if mol is None or mol.GetNumAtoms() == 0:
        return False
    if mol.GetNumConformers() == 0:
        return False
    conf = mol.GetConformer(0)
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        if not (math.isfinite(p.x) and math.isfinite(p.y) and math.isfinite(p.z)):
            return False
    return True

def _compute_gasteiger_ok(mol) -> bool:
    """Compute Gasteiger charges on the given mol (as-is). Return True if all finite."""
    try:
        # Clear any previous charges to avoid stale props
        for a in mol.GetAtoms():
            if a.HasProp('_GasteigerCharge'):
                a.ClearProp('_GasteigerCharge')
            if a.HasProp('_GasteigerH'):
                a.ClearProp('_GasteigerH')
        AllChem.ComputeGasteigerCharges(mol)
    except Exception:
        return False
    for a in mol.GetAtoms():
        if not a.HasProp('_GasteigerCharge'):
            return False
        try:
            q = float(a.GetProp('_GasteigerCharge'))
        except Exception:
            return False
        if not math.isfinite(q):
            return False
    return True

def charges_finite_guard(mol) -> bool:
    """
    Fast check (no hydrogens). If it fails, retry on a copy with explicit H’s.
    We DO NOT write the H-augmented molecule; this is only a screening step.
    """
    if _compute_gasteiger_ok(mol):
        return True
    # retry with hydrogens on a copy
    try:
        molH = Chem.AddHs(Chem.Mol(mol), addCoords=True)
    except Exception:
        return False
    return _compute_gasteiger_ok(molH)

def iter_sdf(input_path, sanitize=True):
    # Forward supplier streams huge files safely
    suppl = Chem.ForwardSDMolSupplier(input_path, sanitize=sanitize, removeHs=False)
    for mol in suppl:
        yield mol

def main():
    ap = argparse.ArgumentParser(
        description="Filter out SDF entries with bad coords or non-finite Gasteiger charges; also drop specific names."
    )
    ap.add_argument("input_sdf")
    ap.add_argument("output_sdf")
    ap.add_argument("--bad-log", default="bad_mols.csv", help="CSV log of removed molecules")
    ap.add_argument("--drop", nargs="*", default=["ZINC000169831752_1"],
                    help="Exact _Name entries to remove (space-separated)")
    ap.add_argument("--no-sanitize", action="store_true",
                    help="Disable RDKit sanitization (faster, but may reduce charge success)")
    ap.add_argument("--progress-every", type=int, default=50000, help="Print progress every N mols")
    args = ap.parse_args()

    w = Chem.SDWriter(args.output_sdf)
    if w is None:
        print("ERROR: Could not open output SDF for writing.", file=sys.stderr)
        sys.exit(1)

    fout = open(args.bad_log, "w", newline="")
    clog = csv.writer(fout)
    clog.writerow(["index_1based", "name", "reason"])

    kept = 0
    dropped = 0
    idx = 0

    drop_set = set(args.drop or [])
    for mol in iter_sdf(args.input_sdf, sanitize=not args.no_sanitize):
        idx += 1
        if idx % max(1, args.progress_every) == 0:
            print(f"...processed {idx:,} molecules (kept {kept:,}, dropped {dropped:,})", file=sys.stderr)

        if mol is None:
            dropped += 1
            clog.writerow([idx, "", "rdkit_parse_failed"])
            continue

        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        # Rule 1: explicit name blocklist
        if name in drop_set:
            dropped += 1
            clog.writerow([idx, name, "name_blocklist"])
            continue

        # Rule 2: finite coordinates present
        if not has_finite_coords(mol):
            dropped += 1
            clog.writerow([idx, name, "bad_or_missing_coords"])
            continue

        # Rule 3: all Gasteiger charges finite (with retry on H-added copy)
        if not charges_finite_guard(mol):
            dropped += 1
            clog.writerow([idx, name, "non_finite_gasteiger"])
            continue

        # Passed all checks → keep
        w.write(mol)
        kept += 1

    w.close()
    fout.close()
    print(f"Done. Processed {idx:,} molecules. Kept {kept:,}, dropped {dropped:,}.")
    print(f"Bad molecules logged to: {args.bad_log}")

if __name__ == "__main__":
    main()
