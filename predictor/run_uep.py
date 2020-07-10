from predictor.make_uep_run import uep

def run_uep(pdb_path, interface_chains, uep_contact_matrix):
    source()
    group1, group2 = check_user_chains(pdb_path, interface_chains)
    uep(pdb_path, uep_contact_matrix, group1, group2)

def check_user_chains(pdb_path, interface_chains):
    group1, group2 = interface_chains.split(",")
    group1_set, group2_set = set([chain for chain in group1]), set([chain for chain in group2])
    pdb_chains = find_chains_in_pdb(pdb_path)
    
    if not group1_set.issubset(pdb_chains) or not group2_set.issubset(pdb_chains):
        raise ValueError("Some stated interface chains are not found in the {} pdb.\n--> Your input: {},{}\n--> Found chains: {}".format(pdb_path, group1, group2, pdb_chains))    
    if len(group1) != len(group1_set) or len(group2) != len(group2_set):
        raise ValueError("Some chains are duplicated within the same group. Make them unique on each group.\n--> Your input: Group1={}, Group2={}".format(group1, group2))
    if bool(group1_set & group2_set):
        raise ValueError("Some chains are repeated in both groups: {}. Do not repeat chains between groups.\n--> Your input: Group1={}, Group2={}".format(group1_set & group2_set, group1, group2))
    print("All stated interface chains {},{} have been found in the pdb {}".format(group1, group2, pdb_path))
    return group1_set, group2_set

def find_chains_in_pdb(pdb_path):
    pdb_chains = set()
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
                pdb_chains.add(chain)
    return pdb_chains

def source():
    print("   -------------------------------   UEP   --------------------------------")
    print("                             Thanks for using UEP")
    print("                      If you used UEP, please, cite us:")
    print("                 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    print("   ------------------------------------------------------------------------\n")
