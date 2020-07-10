import prody
import itertools
from os import path
import numpy as np
prody.confProDy(verbosity="none")

def uep(pdb_path, uep_contact_matrix, group1, group2):
    data = {}
    pdb = prody.parsePDB(pdb_path)
    interface_data = []
    for chain in group1.union(group2):
        interface_contacts = pdb.select("(ca same residue as within 5.00 of (noh chain {0})) and not chain {0}".format(chain))
        for residue_contact, chain_contact in zip(interface_contacts.getResnums().tolist(), interface_contacts.getChids().tolist()):
            if (chain in group1 and chain_contact in group2) or (chain in group2 and chain_contact in group1):
                interface_data.append((residue_contact, chain_contact))
    for interface in interface_data:
        residue, chain = interface[0], interface[1]
        print(residue, chain)
