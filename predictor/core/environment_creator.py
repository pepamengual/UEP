def contacts(folder, radius):
    import glob
    import os
    import prody
    import itertools

    observed_contacts_dictionary = {}
    radius = str(radius)
    pdb_list = glob.glob(os.path.join(folder, "*.pdb"))
    for pdb_file in pdb_list:
        pdb = prody.parsePDB(pdb_file)
        pdb_chains_list = pdb.getChids()

        for chain in set(pdb_chains_list):
            chain_interface = pdb.select('(ca same residue as within '+ radius +'.00 of (noh chain '+ chain +')) and not chain '+ chain +'')
            if chain_interface is None:
                continue

            for residue_id, number_id, chain_id, icode_id, atom_id in zip(chain_interface.getResnames(), chain_interface.getResnums(), chain_interface.getChids(), chain_interface.getIcodes(), chain_interface.getIndices().tolist()):
                if not icode_id:
                    icode_id = "_"
                number_id = '`{}{}`'.format(str(number_id), icode_id)
                first_residue_picker = pdb.select('ca chain '+ chain_id +' and resid '+ number_id +'')
                first_residue_residue_id = tuple(first_residue_picker.getResnames())
                
                residue_interface = pdb.select('(ca same residue as within '+ radius +'.00 of (noh resid '+ number_id +' and chain '+ chain_id +')) and not chain '+ chain_id +'')

                if residue_interface is None:
                    continue

                residue_id_contact_tuple_list = tuple(sorted(residue_interface.getResnames().tolist()))
                residue_id_contact_list = (residue_interface.getResnames().tolist())
                number_id_contact_list = (residue_interface.getResnums().tolist())
                chain_id_contact_list = (residue_interface.getChids().tolist())
                index_id_contact_list = (residue_interface.getIndices().tolist())

                if len(residue_id_contact_tuple_list) < 2:
                    continue

                sorted_contacts = sorted(tuple(zip(residue_id_contact_list, number_id_contact_list, chain_id_contact_list, index_id_contact_list)), key=lambda x:x[0])
                pair_combinatory = (tuple(itertools.combinations(sorted_contacts, 2)))
                residue_id = (first_residue_residue_id)
                for pair in pair_combinatory:
                    combination_contact_residues = (pair[0][0], pair[1][0])
                    if combination_contact_residues in observed_contacts_dictionary:
                        if residue_id in observed_contacts_dictionary[combination_contact_residues]:
                            observed_contacts_dictionary[combination_contact_residues][residue_id] += 1
                        else:
                            observed_contacts_dictionary[combination_contact_residues][residue_id] = 1
                    else:
                        observed_contacts_dictionary[combination_contact_residues] = {residue_id:1}
    return observed_contacts_dictionary
