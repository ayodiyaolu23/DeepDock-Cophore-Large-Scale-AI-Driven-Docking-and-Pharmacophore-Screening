#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
import os
import sys
import time
import gc
from typing import Optional, Tuple

# === CONFIG ===
INPUT_FILE = "2train_smiles_final_updated.smi"   # SMILES per line; optionally "SMILES<tab>NAME"
OUTPUT_DIR = "2train_conformers"                 # Where SDFs will be written
NUM_CONFORMERS = 1                             # Conformers per molecule
NUM_DATABASES = 1                              # Split outputs across N SDF files
FORCE_FIELD = "MMFF"                           # "MMFF" or "UFF"
CHECKPOINT_EVERY = 5000                        # Save progress every N molecules
FLUSH_EVERY = 2000                             # Flush writers & run GC every N molecules
SEED = 17                                      # Set to None for non-deterministic

# === DERIVED PATHS ===
CKPT_PATH = os.path.join(OUTPUT_DIR, ".conformer_ckpt.txt")
LOG_PATH  = os.path.join(OUTPUT_DIR, "conformer_log.txt")


def ensure_output_dir(path: str):
    os.makedirs(path, exist_ok=True)


def load_checkpoint() -> int:
    """
    Returns the next line index to start from (0-based).
    """
    try:
        with open(CKPT_PATH, "r") as f:
            val = f.read().strip()
            return int(val)
    except Exception:
        return 0


def save_checkpoint(idx: int):
    tmp = CKPT_PATH + ".tmp"
    with open(tmp, "w") as f:
        f.write(str(idx))
    os.replace(tmp, CKPT_PATH)


def log(msg: str):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    try:
        with open(LOG_PATH, "a") as f:
            f.write(line + "\n")
    except Exception:
        pass


def parse_smiles_line(line: str) -> Optional[Tuple[str, Optional[str]]]:
    """
    Accepts lines like:
      SMILES
      SMILES<TAB>NAME
      SMILES SPACE NAME  (robust to some whitespace)
    Returns (smiles, name or None). Returns None for blank/invalid lines.
    """
    s = line.strip()
    if not s:
        return None
    # Prefer tab; fallback to whitespace split
    if "\t" in s:
        parts = s.split("\t")
    else:
        parts = s.split()
    if len(parts) == 0:
        return None
    smiles = parts[0]
    name = parts[1] if len(parts) > 1 else None
    return smiles, name


def get_etkdg_params(seed: Optional[int] = None):
    params = AllChem.ETKDGv3()
    if seed is not None:
        params.randomSeed = seed
    # Use all available threads (RDKit interprets 0 as all cores)
    params.numThreads = 0
    # Keep default pruning for first attempt
    return params


def embed_conformers(mol: Chem.Mol, num_confs: int, seed: Optional[int]) -> list:
    """
    Try to embed conformers. If none, retry with more generous params/maxAttempts.
    Returns list of conf IDs (possibly empty).
    """
    params = get_etkdg_params(seed)
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

    if len(conf_ids) == 0:
        # Retry: more attempts + mild pruning to help distinct solutions
        # Handle RDKit version differences for maxAttempts
        try:
            params.SetMaxAttempts(1000)  # Newer RDKit
        except AttributeError:
            pass  # Older builds will ignore; we also pass maxAttempts below
        params.pruneRmsThresh = 0.1
        conf_ids = AllChem.EmbedMultipleConfs(
            mol,
            numConfs=max(1, num_confs),
            params=params,
            maxAttempts=1000  # Works on builds that accept it as a kwarg
        )
    return list(conf_ids)


def optimize_conformers(mol: Chem.Mol, conf_ids: list, force_field: str = "MMFF"):
    if force_field.upper() == "MMFF":
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        if props is None:
            # Fallback to UFF if MMFF params unavailable
            for cid in conf_ids:
                AllChem.UFFOptimizeMolecule(mol, confId=cid)
            return
        for cid in conf_ids:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid, mmffVariant="MMFF94s")
    else:
        for cid in conf_ids:
            AllChem.UFFOptimizeMolecule(mol, confId=cid)


def generate_conformers_for_smiles(smiles: str,
                                   name: Optional[str],
                                   num_confs: int,
                                   force_field: str,
                                   seed: Optional[int]) -> Optional[Chem.Mol]:
    """
    Returns an RDKit molecule with embedded/optimized conformers, or None on failure.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)

    conf_ids = embed_conformers(mol, num_confs, seed)
    if len(conf_ids) == 0:
        return None

    optimize_conformers(mol, conf_ids, force_field=force_field)

    # Annotate props for SDF
    mol.SetProp("_Name", name if name else smiles)
    mol.SetProp("ORIGINAL_SMILES", smiles)
    if name:
        mol.SetProp("TITLE", name)

    return mol


def main():
    ensure_output_dir(OUTPUT_DIR)

    # Prepare writers
    writers = [SDWriter(os.path.join(OUTPUT_DIR, f"database_{i + 1}.sdf"))
               for i in range(NUM_DATABASES)]

    # Load lines
    try:
        with open(INPUT_FILE, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        log(f"ERROR: Input file not found: {INPUT_FILE}")
        sys.exit(1)

    total = len(lines)
    start_idx = load_checkpoint()
    if start_idx > 0:
        log(f"Resuming from checkpoint at index {start_idx} (of {total}).")
    else:
        log(f"Starting fresh. Total molecules: {total}")

    processed_since_flush = 0

    try:
        for i in range(start_idx, total):
            raw = lines[i]
            parsed = parse_smiles_line(raw)
            if parsed is None:
                # skip blank/weird line
                if (i + 1) % 1000 == 0:
                    log(f"Skipped/blank lines up to {i + 1}/{total}")
                continue

            smiles, name = parsed

            try:
                mol = generate_conformers_for_smiles(
                    smiles=smiles,
                    name=name,
                    num_confs=NUM_CONFORMERS,
                    force_field=FORCE_FIELD,
                    seed=SEED
                )
                if mol is None:
                    log(f"FAILED (embed): idx={i} smiles={smiles}")
                else:
                    # Distribute across SDFs by global index to keep deterministic split
                    target = i % NUM_DATABASES
                    # Write each conformer
                    nconf = mol.GetNumConformers()
                    for cid in range(nconf):
                        writers[target].write(mol, confId=cid)
            except Exception as e:
                log(f"ERROR idx={i} smiles={smiles} :: {repr(e)}")

            # Progress & housekeeping
            if (i + 1) % 1000 == 0:
                log(f"Progress: {i + 1}/{total} molecules processed.")

            processed_since_flush += 1
            if processed_since_flush >= FLUSH_EVERY:
                for w in writers:
                    try:
                        w.flush()
                    except Exception:
                        pass
                gc.collect()
                processed_since_flush = 0

            if (i + 1) % CHECKPOINT_EVERY == 0:
                save_checkpoint(i + 1)

        # Final checkpoint to mark completion
        save_checkpoint(total)
        log(f"DONE. Processed {total} molecules. Files saved in '{OUTPUT_DIR}'.")

    finally:
        for w in writers:
            try:
                w.flush()
                w.close()
            except Exception:
                pass


if __name__ == "__main__":
    main()
