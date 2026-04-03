import os

import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
import numpy as np


def convert_cif_to_pdb(cif_path, pdb_path, truncate=True):
    """
    Convert a CIF file to a PDB file using Biotite.
    """
    if not os.path.exists(cif_path):
        print(f"Error: CIF file not found at {cif_path}")
        return

    try:
        # Load the CIF file
        cif_file = pdbx.CIFFile.read(cif_path)

        # Get the structure from the CIF file
        structure = pdbx.get_structure(cif_file, model=1)

        if truncate:
            # PDB standard limits: res_name <= 3, chain_id <= 1
            # Note: Biotite will raise error if these are exceeded
            print(
                "Warning: Truncating long residue names (>3) and chain IDs (>1) for PDB compatibility."
            )

            # Truncate residue names to 3 characters
            new_res_names = np.array([name[:3] for name in structure.res_name])
            structure.res_name = new_res_names

            # Truncate chain IDs to 1 character
            new_chain_ids = np.array([cid[:1] for cid in structure.chain_id])
            structure.chain_id = new_chain_ids

        # Save the structure as a PDB file
        pdb_file = pdb.PDBFile()
        try:
            pdb_file.set_structure(structure)
            pdb_file.write(pdb_path)
            print(f"Successfully converted {cif_path} to {pdb_path}")
        except Exception as inner_e:
            if "exceed" in str(inner_e):
                print(f"Error: PDB format limitation encountered: {inner_e}")
                print(
                    "Hint: PDB format only supports 1-char chain IDs and 3-char residue names."
                )
                print(
                    "Truncation is now ON by default. If you disabled it with --no-truncate, "
                    "re-run without that flag or keep using CIF format."
                )
            else:
                raise inner_e

    except Exception as e:
        print(f"An error occurred during conversion: {e}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert CIF to PDB using Biotite")
    parser.add_argument("input_cif", help="Input CIF file path")
    parser.add_argument("output_pdb", nargs="?", help="Output PDB file path (optional)")
    parser.add_argument(
        "--no-truncate",
        action="store_false",
        dest="truncate",
        help="Do not truncate long residue names and chain IDs",
    )
    parser.set_defaults(truncate=True)

    args = parser.parse_args()

    input_cif = args.input_cif
    output_pdb = (
        args.output_pdb if args.output_pdb else os.path.splitext(input_cif)[0] + ".pdb"
    )

    convert_cif_to_pdb(input_cif, output_pdb, truncate=args.truncate)
