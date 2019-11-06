
def scan_interface(pdb_file, training_data):
    import prody
    residue_list = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    aa_code_triplet_singlet = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V'}

    all_mutant_variants = {}
    pdb = prody.parsePDB(pdb_file)
    pdb_chains = list(set(pdb.getChids()))

    interface_residues_list = []
    for chain in pdb_chains:
        interface_contacts = pdb.select("(ca same residue as within 4.00 of (noh chain {})) and not chain {}".format(chain, chain))
        for residue_contact, chain_contact in zip(interface_contacts.getResnums().tolist(), interface_contacts.getChids().tolist()):
            interface_residues_list.append((residue_contact, chain_contact))

    for contact_tupple in interface_residues_list:
        residue_contact = contact_tupple[0]
        chain_contact = contact_tupple[1]
        residue_selection = pdb.select("ca chain {} and resid {}".format(chain_contact, residue_contact))
        environment_of_residue_selection = pdb.select("(ca same residue as within 4.00 of (noh resid {} and chain {})) and not chain {}".format(residue_contact, chain_contact, chain_contact))
        if environment_of_residue_selection and len(environment_of_residue_selection.getResnums().tolist()) >= 2:
            residue_contact_name = tuple(sorted(residue_selection.getResnames().tolist()))
            environment_of_contact_name = tuple(sorted(environment_of_residue_selection.getResnames().tolist()))
            import itertools
            pair_combination = (tuple(itertools.combinations(environment_of_contact_name, 2)))
            combination_dict = combination_maker(residue_list, pair_combination, training_data, aa_code_triplet_singlet, residue_contact_name, chain_contact, residue_contact)
            for name, score in combination_dict.items():
                all_mutant_variants.setdefault(name, score)
    
    ratio_dict = compute_ratio(all_mutant_variants)
    return ratio_dict

def combination_maker(residue_list, pair_combination, training_data, aa_code_triplet_singlet, residue_contact_name, chain_contact, residue_contact):
    combination_dict = {}
    for target_residue in residue_list:
        target_residue = tuple([target_residue])
        target_score = 0
        for pair in pair_combination:
            pair = tuple(sorted(pair))
            target_score += training_data[pair][target_residue]
            one_letter_original = aa_code_triplet_singlet[residue_contact_name[0]]
            one_letter_target = aa_code_triplet_singlet[target_residue[0]]
            name = "{}{}{}{}".format(one_letter_original, chain_contact, residue_contact, one_letter_target)
            combination_dict.setdefault(name, target_score)
    return combination_dict

def compute_ratio(all_mutant_variants):
    ratio_dict = {}
    for name, score in all_mutant_variants.items():
        original = name[0]
        chain = name[1]
        position = name[2:-1]
        mutation = name[-1]
        original_entry = "{}{}{}{}".format(original, chain, position, original)
        original_score = all_mutant_variants[original_entry]
        ratio = round(score/original_score, 2)
        if ratio > 1:
            original_information = "{}{}{}".format(original, chain, position)
            ratio_dict.setdefault(original_information, {}).setdefault(mutation, ratio)
    return ratio_dict
