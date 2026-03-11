#!/usr/bin/env python3

"""
Author: Yaqi Jiao
Date: 9th March, 2026

mtDNA_tree17_build.py
----------------

Description: 
    This script parses the visually formatted, indentation-based mtDNA phylogenetic tree (Build 17) 
    and converts it into a structured, machine-readable edge list. 
    It performs the following operations:
    1. Reads the raw CSV file where the depth of a node in the tree is represented by the number of preceding empty cells.
    2. Utilizes a Last-In-First-Out (LIFO) Stack logic to track the current path and dynamically determine the direct parent of each new line.
    3. Distinguishes between named haplogroups and unnamed intermediate nodes (branches that only contain mutations).
    4. Constructs a Directed Graph (DiGraph) in memory using NetworkX to map all parent-child relationships and their connecting mutations.
    5. Exports the topological structure (edges) to a flat TSV file for downstream distance calculations.

Input:
    A raw CSV file (e.g., mtDNA_tree_Build_17.csv) where columns represent hierarchical depth (indentation) and cells contain haplogroup names or mutation strings.

Output:
    A tab-separated values file (parsed_mtdna_edges.tsv) detailing the phylogenetic network. 
    Columns include 'Parent_Node', 'Child_Node', and 'Mutations_from_Parent'.

Usage Example:
    $ python mtDNA_tree17_build.py
    # The script will print the total number of nodes constructed and the path to the saved edge list file.
"""

import networkx as nx
import pandas as pd
import re
import csv 
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
input_tree_file = os.path.join(BASE_DIR, "Data", "mtDNA_tree_Build_17.csv")
output_tree_file = os.path.join(BASE_DIR, "Data", "parsed_mtdna_edges.tsv")


print(f"Parsing tree {input_tree_file}...")
tree = nx.DiGraph()  # create nx class, digraph means this is a directed graph, in which the nodes are connected by arrows
tree.add_node("mtDNA_Root", mutations=[])  # set root node

# initialize stack
stack = [("mtDNA_Root", -1)]  # stack is a list, following 'Last-In, First-Out' rule
unnamed_counter = 0 


with open(input_tree_file, 'r', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    
    for row in reader:

        # skip empty row
        if not any(cell.strip() for cell in row): 
            continue
        """
        for each non-empty line in the csv:
        1. Calculate Depth
        2. extract all non-empty cell content, determine mutation
        3. Find Parent (Stack operation):
            Continuously compare the depth of the new row with the last node in the stack.
            * WHILE stack_top_depth >= current_depth: pop (discard) the last node from the stack. 
            * The node that remains at the top of the stack is guaranteed to be the direct Parent.
        4. Update Graph: Add the new node to the Tree/Graph, and draw an edge between the Parent and this new node.
        5. Append the new node's name and depth into the stack for future rows to reference.
        """
        # calculate depth thorugh counting empty cells
        depth = 0
        for cell in row:
            if cell.strip() == "":
                depth += 1
            else:
                break # Stop counting when the first cell containing a word (haplogroup or mutation) is encountered

        # A token is the smallest meaningful block of text extracted from a long string of meaningless characters    
        tokens = [cell.strip() for cell in row if cell.strip() != ""]  # collect all non-empty cells
        if not tokens: continue
        
        first_token = tokens[0]
        
        # Determine whether the first word is a haplogroup name or a mutation.
        is_mutation = bool(re.match(r'^\(?[ACGT]\d+', first_token.upper()))
        
        if is_mutation:
            hg_name = f"Unnamed_Branch_{unnamed_counter}"
            unnamed_counter += 1
            mutations = tokens 
        else:
            hg_name = first_token
            mutations = tokens[1:] 
            
        # Find the direct parent node (core pop/push logic)
        while stack and stack[-1][1] >= depth:
            stack.pop()
            
        parent_name = stack[-1][0]
        
        # Add the current node to the graph and connect it to its parent node.
        tree.add_node(hg_name, mutations=mutations)
        tree.add_edge(parent_name, hg_name)
        
        # Put onto the stack
        stack.append((hg_name, depth))

print(f"Tree constructed! {tree.number_of_nodes()} nodes in total")


print("Saving tree structure...")
edge_data = []

for parent, child in tree.edges():
    muts = tree.nodes[child].get('mutations', [])
    # use ',' to connect mutation
    muts_str = ",".join(muts) if muts else ""
    
    edge_data.append({
        "Parent_Node": parent,
        "Child_Node": child,
        "Mutations_from_Parent": muts_str
    })

df_edges = pd.DataFrame(edge_data)
df_edges.to_csv(output_tree_file, sep='\t', index=False)

print(f"Parsing finished! Tree file saved in: {output_tree_file}")