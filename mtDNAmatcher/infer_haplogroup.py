#!/usr/bin/env python3

"""
Author: Yaqi Jiao
Date: 12th March, 2026

infer_haplogroup.py
----------------
Description:
    This module provides the core inference engine for determining a user's maternal 
    haplogroup from raw mtDNA variants. It employs a bottom-up cumulative scoring 
    algorithm against the PhyloTree (Build 17) topology.

    Key functionalities:
    1. Loads the parsed phylogenetic tree into an optimized dictionary structure.
    2. Dynamically accumulates all ancestral mutations for any given haplogroup 
    by traversing upwards to the mtDNA Root.
    3. Handles complex PhyloTree mutation nomenclatures (e.g., back mutations '!', 
    uncertain mutations '()') via regular expression parsing.
    4. Calculates a robust similarity score between the user's SNPs and the haplogroup's 
    defining SNPs. Critically, the scoring mechanism is resistant to missing data, 
    as it normalizes the score only against the positions actively tested in the user's file.
"""


import pandas as pd
import re

def load_mutation_tree_from_tsv(tsv_path):
    """
    1. Read the pre-parsed TSV file into a pandas DataFrame.
    2. Iterate through each row to extract 'Child_Node', 'Parent_Node', and 'Mutations_from_Parent'.
    3. Clean the mutation string: if it's NaN or empty, assign an empty list
    4. Use re.split to split the mutation string into a list
    5. Map the child node as the dictionary key, and store a tuple of (parent_node, mutation_list) as its value.
    """
    tree_dict = {}
    df = pd.read_csv(tsv_path, sep='\t')
    for _, row in df.iterrows():
        child = str(row['Child_Node'])
        parent = str(row['Parent_Node'])
        # compatible with comma or space seprated mutation information
        mut_str = str(row['Mutations_from_Parent'])
        if pd.isna(row['Mutations_from_Parent']) or mut_str.strip() == "":
            muts = []
        else:
            muts = re.split(r'[,\s]+', mut_str.strip())
        tree_dict[child] = (parent, muts)
    return tree_dict

# split the data into position and genotype, and remove non-standard values
def parse_tree_mutation(mut_str):
    clean_str = mut_str.strip('!()')

    match = re.match(r'([a-zA-Z]+)(\d+)([a-zA-Z]+)', clean_str)
    if match:
        return int(match.group(2)), match.group(3).upper()
    return None, None

def get_haplogroup_snps(target_haplogroup, tree_dict):
    all_snps = set()  # Initialize an empty Set (`all_snps`) to store unique mutations and prevent duplicates
    current_node = target_haplogroup  # Set the starting point (`current_node`) to the user's target haplogroup

    # Enter a While loop: as long as the current node exists in the dictionary, fetch its parent and mutations
    while current_node in tree_dict:
        parent, muts = tree_dict[current_node]
        all_snps.update(muts)
        current_node = parent  # Update `current_node` to be the parent node
        if current_node == 'mtDNA_Root': break
    return all_snps

def infer_best_haplogroup(user_snps, tree_dict):
    print("Calculating similarities across the phylogenetic tree...")
    results = []
    """
    Calculate the similarity between the user's DNA and each haplogroup on the tree
    1. For each haplogroup in the tree, retrieve its complete list of accumulated mutations from the Root
    2. if there is no mutation in this line, try next haplogroup
    3. if there are mutation information: Parse these defining mutations into specific genomic positions and derived alleles
    4. Compare against the user's provided SNPs, count tested position and matched position numbers
    5. If the position was tested by the user,check if the derived allele matches
    6. Calculate a confidence score (matched / valid tested positions), store them in a dictionary
    7. return the Top 5 matches
    """
    # Iterate every haplogroup on the tree
    for haplo in tree_dict.keys():

        defining_snps = get_haplogroup_snps(haplo, tree_dict)
        if not defining_snps: continue
            
        tested_positions = 0  
        matched_positions = 0 
        
        for mut in defining_snps:
            pos, derived_alt = parse_tree_mutation(mut)
            if pos is None: continue
            
            if pos in user_snps:
                tested_positions += 1
                if user_snps[pos] == derived_alt:
                    matched_positions += 1
        
        if tested_positions > 0:
            score = matched_positions / tested_positions
            results.append({
                'haplogroup': haplo,
                'score': score,
                'matched': matched_positions,
                'tested': tested_positions
            })
            
    # Sort by score; if scores are the same, choose the one with more absolute matches.
    results.sort(key=lambda x: (x['score'], x['matched']), reverse=True)
    return results[:5]

# if __name__ == '__main__':
#     BASE_DIR = os.path.dirname(os.path.abspath(__file__))
#     tree_tsv_file = os.path.join(BASE_DIR, "Data", "parsed_mtdna_edges.tsv")
#     user_test_file = os.path.join(BASE_DIR, "Data", "23andme_v2_test.txt")

#     print("\n" + "="*50)
#     print("mtDNA Matcher Pipeline Started")
#     print("="*50 + "\n")

#     print(">>> STEP 1: Loading Phylogenetic Tree...")
#     my_tree_dict = load_mutation_tree_from_tsv(tree_tsv_file)
#     print(f"Loaded {len(my_tree_dict)} haplogroups from the database.\n")

#     print(">>> STEP 2: Processing User DNA Data...")
#     if not os.path.exists(user_test_file):
#         print(f"Error: Cannot find user data at {user_test_file}")
#         print("Please make sure you saved the test data there!")
#         exit()
#     my_user_snps = parse_user_dna(user_test_file)
#     print("")

#     print(">>> STEP 3: Inferring Maternal Haplogroup...")
#     top_matches = infer_best_haplogroup(my_user_snps, my_tree_dict)
    
#     print("\n --- Inference Results --- ")
#     for i, match in enumerate(top_matches):
#         print(f" Top {i+1}: Haplogroup {match['haplogroup']}")
#         print(f" Confidence: {match['score']:.1%} ({match['matched']} out of {match['tested']} tested SNPs matched)\n")

#     best_haplogroup = top_matches[0]['haplogroup']
#     print(f"Conclusion: The user's haplogroup is automatically inferred as: **{best_haplogroup}**")
#     print(">>> Next step: Route this haplogroup to the Ancient Database (Your previous code)!")
#     print("\n" + "="*50 + "\n")