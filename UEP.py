import argparse
import os
from compress_pickle import load
from predictor.uep_main import run_uep

HELP = " \
Command:\n \
----------\n \
UEP scans the protein-protein interface of a PDB file: \n \
    python3 UEP.py --pdb=pdb_path.pdb interface=Group1,Group2\
"

def parse_args():
    """
    """
    parser = argparse.ArgumentParser(description=HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pdb", type=str, help="Path to the pdb to be evaluated", default="")
    parser.add_argument("--interface", type=str, help="'Group1,Group2' --> 'AB,C' or 'ABC,DEFG'", default="")
    args = parser.parse_args()
    return args.pdb, args.interface

def main(pdb_path="", interface_chains=""):
    """
    """
    if not pdb_path or not os.path.exists(pdb_path):
        raise ValueError("--> Please, define a valid pdb path to score")
    if not interface_chains or len(interface_chains.split(",")) != 2:
        raise ValueError("--> Please, define an interface chain following the scheme: Group1,Group2")
    
    uep_contact_matrix = load("interactome_model/model.lzma", compression="lzma", set_default_extension=False)   
    uep_predictions = run_uep(pdb_path, interface_chains, uep_contact_matrix)

if __name__ == "__main__":
    pdb, interface = parse_args()
    main(pdb_path=pdb, interface_chains=interface)
