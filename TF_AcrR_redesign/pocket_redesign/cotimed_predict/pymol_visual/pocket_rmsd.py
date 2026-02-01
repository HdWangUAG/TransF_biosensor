#!/usr/bin/env python3
import argparse
import sys

# Use PyMOL in quiet, no-GUI mode if run as "python script.py"
try:
    import pymol
    pymol.finish_launching(["pymol", "-qc"])
except Exception:
    pass

from pymol import cmd

def parse_resi_spec(spec):
    """
    Parse pocket residue spec like:
      "A:45,46,82,110-115;B:10-12,30"
    Returns dict: {(chain) -> list of resi strings}
    """
    per_chain = {}
    if not spec:
        return per_chain
    # Split chain blocks by ';'
    blocks = [b.strip() for b in spec.split(";") if b.strip()]
    for block in blocks:
        if ":" not in block:
            raise ValueError(f"Pocket residue block must be 'CHAIN:resi_list', got '{block}'")
        chain, resi_list = block.split(":", 1)
        chain = chain.strip()
        # Accept commas or + separators
        items = [r.strip() for r in resi_list.replace("+", ",").split(",") if r.strip()]
        per_chain.setdefault(chain, []).extend(items)
    return per_chain

def build_resi_selection_from_map(obj, chain_to_resi):
    """
    Build a PyMOL selection string for residues defined by chain_to_resi on object 'obj'.
    """
    parts = []
    for chain, resis in chain_to_resi.items():
        parts.append(f"({obj} and chain {chain} and resi " + "+".join(resis) + ")")
    if not parts:
        return f"{obj} and none"
    return "(" + " or ".join(parts) + ")"

def residue_keys_from_selection(sel):
    """
    Return set of (chain, resi) for residues present in selection 'sel'.
    """
    model = cmd.get_model(sel)
    keys = set()
    for a in model.atom:
        keys.add((a.chain, a.resi))
    return keys

def build_selection_from_residue_keys(obj, keys):
    """
    Build selection on 'obj' from a set of (chain, resi) tuples.
    """
    per_chain = {}
    for chain, resi in keys:
        per_chain.setdefault(chain, []).append(resi)
    return build_resi_selection_from_map(obj, per_chain)

def atom_mask_from_choice(choice):
    """
    Map atom set choice to a PyMOL atom mask applied within the pocket selection.
    """
    if choice == "heavy":
        return "not name H*"
    if choice == "backbone":
        return "name N+CA+C+O"
    if choice == "ca":
        return "name CA"
    if choice == "all":
        return "all"
    raise ValueError(f"Unknown atom-set: {choice}")

def main():
    ap = argparse.ArgumentParser(description="Compute pocket-only RMSD between two structures using PyMOL.")
    ap.add_argument("--ref", required=True, help="Reference structure (PDB/mmCIF/etc.)")
    ap.add_argument("--mob", required=True, help="Mobile structure to compare")
    ap.add_argument("--object-names", nargs=2, default=["ref", "mob"], help="Names to use for loaded objects")
    ap.add_argument("--pocket-resi", default=None,
                    help="Pocket by residue list: e.g., 'A:45,46,82,110-115;B:10-12'")
    ap.add_argument("--lig-sel", default=None,
                    help="Alternative pocket: residues within radius of this ligand selection, e.g., 'resn LIG and chain X'")
    ap.add_argument("--radius", type=float, default=5.0, help="Radius (Å) for ligand-based pocket (default: 5.0)")
    ap.add_argument("--atom-set", choices=["heavy", "backbone", "ca", "all"], default="backbone",
                    help="Atoms to include in RMSD (default: backbone)")
    ap.add_argument("--no-global-align", action="store_true", help="Skip global Cα alignment before measuring")
    ap.add_argument("--save-aligned", default=None, help="Optional: save aligned mobile as PDB")
    args = ap.parse_args()

    ref_name, mob_name = args.object_names

    # Load structures
    cmd.reinitialize()
    cmd.load(args.ref, ref_name)
    cmd.load(args.mob, mob_name)

    # remove hydrogen on the residues to avoid discrepancies

    # Optional global alignment on Cα atoms
    if not args.no_global_align:
        cmd.align(f"{mob_name} and polymer and name CA", f"{ref_name} and polymer and name CA")

    # Determine pocket selections on ref and mob
    if args.pocket_resi:
        # Pocket defined explicitly by chain:resi list(s)
        chain_map = parse_resi_spec(args.pocket_resi)
        pocket_ref = f"{ref_name} and polymer and byres ({build_resi_selection_from_map(ref_name, chain_map)})"
        pocket_mob = f"{mob_name} and polymer and byres ({build_resi_selection_from_map(mob_name, chain_map)})"
    elif args.lig_sel:
        # Pocket defined by residues within radius of ligand — derive residue set from ref,
        # then map the same residues (by chain+resi) onto mob for consistent pairing.
        lig_ref = f"{ref_name} and ({args.lig_sel})"
        pocket_atoms_ref = f"{ref_name} and polymer within {args.radius} of ({lig_ref})"
        pocket_ref = f"{ref_name} and polymer and byres ({pocket_atoms_ref})"

        # Extract residue IDs from ref pocket and build corresponding mob selection
        keys = residue_keys_from_selection(pocket_ref)
        if not keys:
            print("Error: reference pocket selection is empty. Check --lig-sel and --radius.", file=sys.stderr)
            sys.exit(2)
        pocket_mob = f"{mob_name} and polymer and byres ({build_selection_from_residue_keys(mob_name, keys)})"
    else:
        print("Error: specify either --pocket-resi or (--lig-sel and --radius).", file=sys.stderr)
        sys.exit(2)

    # Apply atom mask inside pocket
    atom_mask = atom_mask_from_choice(args.atom_set)
    sel_ref = f"{pocket_ref} and ({atom_mask})"
    sel_mob = f"{pocket_mob} and ({atom_mask})"

    # Sanity checks
    n_ref = cmd.count_atoms(sel_ref)
    n_mob = cmd.count_atoms(sel_mob)
    if n_ref == 0 or n_mob == 0:
        print(f"Error: empty selection. ref atoms={n_ref}, mob atoms={n_mob}.", file=sys.stderr)
        sys.exit(2)
    if n_ref != n_mob:
        # Often harmless if atoms map consistently, but warn user.
        print(f"Warning: atom counts differ (ref={n_ref}, mob={n_mob}). RMSD may be unreliable.", file=sys.stderr)

    # Compute RMSD without further fitting (use result after global alignment)
    rmsd = cmd.rms_cur(sel_mob, sel_ref)
    print(f"Pocket RMSD ({args.atom_set}) = {rmsd:.3f} Å")

    # Optionally save aligned mobile
    if args.save_aligned:
        cmd.save(args.save_aligned, mob_name)
        print(f"Saved aligned mobile to: {args.save_aligned}")

if __name__ == "__main__":
    main()
