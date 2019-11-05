def score_skempi(training_data, skempi_pdb_file, radius):
    import prody
    import itertools

    radius = str(radius)
    pdb_file_name, pdb_code, resume_mutation, original_residue, chain, position, mutant_residue = define_mutation(skempi_pdb_file)
    original_residue_triplet, mutant_residue_triplet = singlets_to_triplets(original_residue, mutant_residue)

    pdb = prody.parsePDB(skempi_pdb_file)
    environment_residue = pdb.select('(ca same residue as within '+ radius +'.00 of (noh resid '+ position +' and chain '+ chain +')) and not chain '+ chain +'')
    if environment_residue is None:
        print("Position {} of {} pdb file has no contacts closer than {}.00 radius.".format(position, pdb_file_name, radius))
        return
    select_original_residue = pdb.select('ca chain '+ chain +' and resid '+ position +'') # Should be the same than original_residue_triplet
    target_original_residue = tuple(sorted(select_original_residue.getResnames().tolist()))
    target_mutant_residue = (mutant_residue_triplet,)
    residue_id_contact_tuple_list = tuple(sorted(environment_residue.getResnames().tolist()))

    if len(residue_id_contact_tuple_list) < 2:
        print("Position {} of {} pdb file has less than 2 contacts and can not be scored.".format(position, pdb_file_name))
        return

    pair_combinatory = (tuple(itertools.combinations(residue_id_contact_tuple_list, 2)))
    
    contact_count_original_residue = 0
    contact_count_mutant_residue = 0

    for pair in pair_combinatory:
        environment_of_mutation = tuple(sorted(pair))
        if environment_of_mutation in training_data:
            contact_count_original_residue += training_data[environment_of_mutation][target_original_residue]
            contact_count_mutant_residue += training_data[environment_of_mutation][target_mutant_residue]
    name = "{}_{}".format(pdb_code, resume_mutation)
    ratio = round(contact_count_mutant_residue/contact_count_original_residue)
    return name, ratio

def define_mutation(skempi_pdb_file):
    pdb_file_name = skempi_pdb_file.split("/")[-1]
    pdb_code = pdb_file_name.split(".pdb")[0].split("_")[1]
    resume_mutation = pdb_file_name.split(".pdb")[0].split("_")[-1].split("-")[0]
    original_residue = resume_mutation[0]
    chain = resume_mutation[1]
    position = resume_mutation[2:-1]
    mutant_residue = resume_mutation[-1]
    return pdb_file_name, pdb_code, resume_mutation, original_residue, chain, position, mutant_residue

def singlets_to_triplets(original_residue, mutant_residue):
    aa_code_singlet_to_triplet = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN','G':'GLY','H':'HIS',
            'I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP',
            'Y':'TYR','V':'VAL'}
    original_residue_triplet = aa_code_singlet_to_triplet[original_residue]
    mutant_residue_triplet = aa_code_singlet_to_triplet[mutant_residue]
    return original_residue_triplet, mutant_residue_triplet
