import prody
import itertools
import pandas as pd
from os import path
import numpy as np
prody.confProDy(verbosity="none")

def uep(pdb_path, uep_contact_matrix, group1, group2):
    aa_code = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS',
               'E':'GLU','Q':'GLN','G':'GLY','H':'HIS','I':'ILE',
               'L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO',
               'S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
    aa_list = list(aa_code.values())
    data = {}
    pdb = prody.parsePDB(pdb_path)
    interface_data = []
    for chain in group1.union(group2):
        interface_contacts = pdb.select("(ca same residue as within 5.00 of (noh chain {0})) and not chain {0}".format(chain))
        for position_contact, resname_contact, chain_contact in zip(interface_contacts.getResnums().tolist(), interface_contacts.getResnames().tolist(), interface_contacts.getChids().tolist()):
            if (chain in group1 and chain_contact in group2) or (chain in group2 and chain_contact in group1):
                interface_data.append((chain_contact, position_contact, resname_contact))

    for interface in sorted(interface_data):
        chain, position, original_residue = interface[0], interface[1], interface[2]
        near_residues_original = pdb.select('(ca same residue as within 5.00 of (noh resid '+ str(position) +' and chain '+ chain +')) and not chain '+ chain +'')
        near_residues_mutation_4, near_residues_mutation_5, near_residues_mutation_6 = approximate_volume(pdb, position, chain)
        for mutant_residue in aa_list:
            if original_residue != mutant_residue:
                counts_original, counts_mutation = 0, 0
                near_residues_mutation = volume(original_residue, mutant_residue, near_residues_mutation_4, near_residues_mutation_5, near_residues_mutation_6)
                if near_residues_original == None or near_residues_mutation == None:
                    continue               
                original, mutant = (original_residue,), (mutant_residue,)
                environment_original = tuple(sorted(near_residues_original.getResnames().tolist()))
                environment_mutation = tuple(sorted(near_residues_mutation.getResnames().tolist()))

                if len(environment_original) < 2 or len(environment_mutation) < 2:
                    continue

                pair_combinatory_original = (tuple(itertools.combinations(environment_original, 2)))
                pair_combinatory_mutation = (tuple(itertools.combinations(environment_mutation, 2)))

                for pair in pair_combinatory_original:
                    environment_of_mutation = tuple(sorted(pair))
                    if environment_of_mutation in uep_contact_matrix and original in uep_contact_matrix[environment_of_mutation]:
                        counts_original += uep_contact_matrix[environment_of_mutation][original] * sum(uep_contact_matrix[environment_of_mutation].values())
                for pair in pair_combinatory_mutation:
                    environment_of_mutation = tuple(sorted(pair))
                    if environment_of_mutation in uep_contact_matrix and mutant in uep_contact_matrix[environment_of_mutation]:
                        counts_mutation += uep_contact_matrix[environment_of_mutation][mutant] * sum(uep_contact_matrix[environment_of_mutation].values())
                if counts_mutation != 0 and counts_original != 0:
                    ratio = (counts_mutation/counts_original)
                    ddG = round(-np.log(ratio), 3)
                    data.setdefault("{0}_{1}_{2}".format(chain, position, original_residue), {}).setdefault(mutant_residue, ddG)
    df = generate_dataframe(data)
    return df

def generate_dataframe(data):
    df = pd.DataFrame.from_dict(data, orient="index")
    df1 = df.index.str.split('_', expand=True).to_frame()
    df1.columns = list('abc')
    df1.index = df.index
    df1[['b']] = df1[['b']].astype(int)
    df1 = df1.sort_values(['a', 'b'])
    df = df.reindex(df1.index)
    return df

def approximate_volume(pdb, position, chain):
    near_residues_mutation_4 = pdb.select('(ca same residue as within 4.00 of (noh resid '+ str(position) +' and chain '+ chain +')) and not chain '+ chain +'')
    near_residues_mutation_5 = pdb.select('(ca same residue as within 5.00 of (noh resid '+ str(position) +' and chain '+ chain +')) and not chain '+ chain +'')
    near_residues_mutation_6 = pdb.select('(ca same residue as within 6.00 of (noh resid '+ str(position) +' and chain '+ chain +')) and not chain '+ chain +'')
    return near_residues_mutation_4, near_residues_mutation_5, near_residues_mutation_6

def volume(aa_original, aa_mutation, near_residues_mutation_4, near_residues_mutation_5, near_residues_mutation_6):
    volume_dict = {"ALA": 88.6, "ARG": 173.4, "ASN": 114.1, "ASP": 111.1, "CYS": 108.5, 
                   "GLN": 143.8, "GLU": 138.4, "GLY": 60.1, "HIS": 153.2, "ILE": 166.7, 
                   "LEU": 166.7, "LYS": 168.6, "MET": 162.9, "PHE": 189.9, "PRO": 112.7, 
                   "SER": 89.0, "THR": 116.1, "TRP": 227.8, "TYR": 193.6, "VAL": 140}
    mutation, original = volume_dict[aa_mutation], volume_dict[aa_original]
    ratio_volume = mutation/original
    if ratio_volume > 1.3:
        return near_residues_mutation_6
    if ratio_volume < 0.7:
        return near_residues_mutation_4
    if ratio_volume > 0.7 and ratio_volume < 1.3:
        return near_residues_mutation_5
