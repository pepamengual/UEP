def training_with_multiprocessing(radius, number_of_processors, path_training_folders):
    import glob
    from multiprocessing import Pool
    from predictor.core import environment_creator

    training_data = {}

    pool = Pool(processes=number_of_processors)
    multiple_results = []
    for folder in glob.glob(path_training_folders):
        multiple_results.append(pool.apply_async(environment_creator.contacts, (folder, radius)))
    for result in multiple_results:
        observed_contacts_dictionary = result.get()
        for pair_residues, target_residue_dictionary in observed_contacts_dictionary.items():
            if pair_residues in training_data:
                for target_residue, observed_contacts in target_residue_dictionary.items():
                    if target_residue in training_data[pair_residues]:
                        training_data[pair_residues][target_residue] += observed_contacts
                    else:
                        training_data[pair_residues][target_residue] = observed_contacts
            else:
                 training_data[pair_residues] = target_residue_dictionary
    pool.terminate()

    return training_data
