import prody
import itertools
from os import path
import numpy as np


def uep(pdb_path, uep_contact_matrix, group1, group2):
    data = {}
    pdb = prody.parsePDB(path_pdb)
    all_chains = group1.union(group2)
    interface_data = []
    for chain in all_chains:
        interface_contacts = pdb.select("(ca same residue as within 5.00 of (noh chain {0})) and not chain {0}".format(chain))
        for residue_contact, chain_contact in zip(interface_contacts.getResnums().tolist(), interface_contacts.getChids().tolist()):
            if (chain in group1 and chain_contact in group2) or (chain in group2 and chain_contact in group2):
                interface_data.append((residue_contact, chain_contact))
    for interface in interface_data:
        residue, chain = interface[0], interface[1]
        print(residue, chain)








def scoring_skempi(pdb_name, candidate_list, training_data):
    aa_code_singlet_to_triplet = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN','G':'GLY','H':'HIS',
            'I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP',
            'Y':'TYR','V':'VAL'}
    frequency_random_model = get_frequency_random_model(training_data)
    data = {}
    radius = str(5.00)
    path_pdb = "PDBs/{}.pdb".format(pdb_name)
    if not path.exists(path_pdb):
        return None
    pdb = prody.parsePDB(path_pdb)
    for candidates in candidate_list:
        candidate_chunks = candidates.split("_")[-1].split(",")
        counts_original, counts_mutation = 0, 0
        list_to_check = {}
        for candidate in candidate_chunks:
            chain = candidate[1]
            aa_original = candidate[0]
            position = candidate[2:-1]
            aa_mutation = candidate[-1]
            list_to_check.setdefault("{}_{}_{}".format(chain, aa_code_singlet_to_triplet[aa_original], position), aa_code_singlet_to_triplet[aa_mutation])
        for candidate in candidate_chunks:
            chain = candidate[1]
            aa_original = candidate[0]
            position = candidate[2:-1]
            aa_mutation = candidate[-1]
            new_ratio = volume_corr(aa_original, aa_mutation)
            #add changing near residues function for all more than one mutations
            near_residues_original = pdb.select('(ca same residue as within '+ radius +' of (noh resid '+ position +' and chain '+ chain +')) and not chain '+ chain +'')
            near_residues_mutation = pdb.select('(ca same residue as within '+ new_ratio +' of (noh resid '+ position +' and chain '+ chain +')) and not chain '+ chain +'')
            if near_residues_original == None or near_residues_mutation == None:
                 continue

            aa_original_three = (aa_code_singlet_to_triplet[aa_original],)
            aa_mutation_three = (aa_code_singlet_to_triplet[aa_mutation],)
        
            environment_original = tuple(sorted(near_residues_original.getResnames().tolist()))
            if len(candidate_chunks) > 1:
                resnames = near_residues_original.getResnames().tolist()
                chains = near_residues_original.getChids().tolist()
                resnums = near_residues_original.getResnums().tolist()
                for i, r in enumerate(resnames):
                    find = "{}_{}_{}".format(chains[i], resnames[i], resnums[i])
                    if find in list_to_check.keys():
                        resnames[i] = list_to_check[find]

                environment_mutation = tuple(sorted(resnames))
            else:
                environment_mutation = tuple(sorted(near_residues_mutation.getResnames().tolist()))
            if len(environment_original) < 2 or len(environment_mutation) < 2: # or
                continue
            pair_combinatory_original = (tuple(itertools.combinations(environment_original, 2)))
            pair_combinatory_mutation = (tuple(itertools.combinations(environment_mutation, 2)))

            for pair in pair_combinatory_original:
                environment_of_mutation = tuple(sorted(pair))
                if environment_of_mutation in training_data and aa_original_three in training_data[environment_of_mutation]:
                    counts_original += training_data[environment_of_mutation][aa_original_three] * sum(training_data[environment_of_mutation].values())
            for pair in pair_combinatory_mutation:
                environment_of_mutation = tuple(sorted(pair))
                if environment_of_mutation in training_data and aa_mutation_three in training_data[environment_of_mutation]:
                    counts_mutation += training_data[environment_of_mutation][aa_mutation_three] * sum(training_data[environment_of_mutation].values())
        if counts_mutation != 0 and counts_original != 0:   
            ratio = (counts_mutation/counts_original)# * (frequency_random_model[aa_mutation_three]/frequency_random_model[aa_original_three])
            #if len(candidates.split("_")[1]) == 1 and len(candidates.split("_")[2]) == 1:
            #print(candidates, ratio)
            ddG = round(-np.log(ratio), 3)
            data.setdefault(candidates, ddG)
    return data

def volume_corr(aa_original, aa_mutation):
    volume_dict = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8, "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140}
    mutation, original = volume_dict[aa_mutation], volume_dict[aa_original]
    min_volume, max_volume = min(volume_dict.values()), max(volume_dict.values())
    diff_volume_min_max = np.abs(max_volume - min_volume)
    diff_volume_mutation_original = np.abs(mutation - original)

    ratio_volume = mutation/original
    if ratio_volume > 1.3:
        new_ratio = str(6)
    if ratio_volume < 0.7:
        new_ratio = str(4)
    if ratio_volume > 0.7 and ratio_volume < 1.3:
        new_ratio = str(5)
    return new_ratio


def get_frequency_random_model(training_data):
    frequency_random_model = {}
    for environment, target_dict in training_data.items():
        for target, count in target_dict.items():
            frequency_random_model.setdefault(target, 0)
            frequency_random_model[target] += count
    return frequency_random_model

def run_multiprocessing(experimental_skempi_ratios, cpus, training_data):
    skempi_uep_predictions = {}
    pool = Pool(processes=cpus)
    multiple_results = []
    
    experimental_skempi_data = {}
    for candidate, aff_ratio in experimental_skempi_ratios.items():
        pdb = candidate.split("_")[0]
        experimental_skempi_data.setdefault(pdb, {}).setdefault(candidate, aff_ratio)

    for pdb, candidate_dict in experimental_skempi_data.items():
        candidate_list = list(candidate_dict.keys())
        multiple_results.append(pool.apply_async(scoring_skempi, (pdb, candidate_list, training_data)))
    for result in multiple_results:
        data = result.get()
        if data != None:
            for name, prediction in data.items():
                skempi_uep_predictions.setdefault(name, prediction)
    pool.terminate()
    return skempi_uep_predictions

