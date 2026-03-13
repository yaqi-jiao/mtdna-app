#!/usr/bin/env python3

"""
Author: Yaqi Jiao
Date: 9th March, 2026

mtDNAmatcher.py
----------------

This application is an end-to-end bioinformatics pipeline that infers a user's maternal haplogroup 
    from raw DNA sequencing data and identifies their closest ancient human relatives.
    
    1. Load background data: 
        - Reconstructs the mtDNA phylogenetic tree (Build 17) as an undirected graph for topological distance calculations.
        - Loads a structured mutation dictionary for rapid haplogroup inference.
        - Loads the cleaned ancient DNA database (AADR).
        
    2. Read user input (Dual Mode): Acquires input from the GUI text entry.
        - Mode A (Raw DNA File): If a file path is provided, the script parses raw genotype data (e.g., 23andMe format), 
        extracts valid mtDNA SNPs, and automatically infers the most likely haplogroup by scoring user mutations against 
        the definitive nodes in the phylogenetic tree.
        - Mode B (Direct Input): If a standard string is provided, it treats the input directly as the target haplogroup.
        
    3. Calculate and Match: 
        Validates the (inferred or provided) haplogroup, calculates the genetic distance (shortest path in the tree graph) 
        between the user and ancient individuals, and sorts them to find the closest evolutionary matches.

Input:
    Via GUI text entry. Accepts either:
    - A string representing a valid mtDNA haplogroup (e.g., "H1a", "X2c2").
    - A valid file path to a raw DNA text file containing genomic coordinates and genotypes (e.g., "TestData/user_test.txt").

Output:
    Displays results in the GUI text area:
    - (If file input): The haplogroup inference process, displaying top candidate haplogroups with confidence scores.
    - The Top 5 closest ancient matches, including their Sample ID, Genetic distance (steps in the phylogenetic tree), 
    Haplogroup, Period of life (Age BP), and Location of origin.

Usage Example:
    $ python mtDNAmatcher.py
    # A GUI window will launch. 
    # Enter "X2" OR a file path like "C:\\TestData\\user_test.txt" into the input box and click "Start Matching...".
"""

import pandas as pd
import networkx as nx
import tkinter as tk
from tkinter import messagebox
import os

from mtDNAmatcher.infer_haplogroup import load_mutation_tree_from_tsv, infer_best_haplogroup
from mtDNAmatcher.user_data_process import parse_user_dna


# Load background data
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
print("Loading data...")
aadr_file = os.path.join(BASE_DIR, "Data", "Clean_metadata_DB.tsv") 
tree_file = os.path.join(BASE_DIR, "Data", "parsed_mtdna_edges.tsv")
tree_df = pd.read_csv(tree_file, sep='\t')

# create empty graph, then iterate tree_df row by row , 
"""
Architecture Note:
In `mtDNA_tree17_build.py`, the phylogenetic tree is built as a Directed Graph (DiGraph) in memory, 
but only the edge pairs are exported to the TSV file (no serialized graph object is saved).

In this script, we reconstruct the tree topology from the TSV file. 
Crucially, we instantiate it as an UNDIRECTED graph (`nx.Graph`) because calculating genetic distance 
requires bidirectional traversal: tracing upwards to the Most Recent Common Ancestor (MRCA) 
and then downwards to a collateral relative.
"""
G = nx.Graph() 
for _, row in tree_df.iterrows():
    G.add_edge(row['Parent_Node'], row['Child_Node'])
valid_hgs = set(G.nodes())
tree_dict = load_mutation_tree_from_tsv(tree_file)

aadr_df = pd.read_csv(aadr_file, sep='\t')
aadr_df = aadr_df[aadr_df['Haplogroup'].isin(valid_hgs)].copy()
print("Data loaded, launch app")

# defination function
def run_match():
    # acquire user input
    input_text = entry_hg.get().strip().strip('"').strip("'")  # Extract the text from the input box on the interface and remove the leading and trailing spaces
    text_result.delete(1.0, tk.END) 
    
    if not input_text:
        messagebox.showwarning("WARNING", "input file path or haplogroup!")
        return

    user_hg = ""  # create empty string to store haplogroup 

    """
    To acquire haplogroup from user input:
    1. Determine if the input is an existing file path
        if the file exists, use [[][]]
    2. If the input is not a file path, then use input as hg_name
    3. Search valid_hgs with user_hg, to make sure the haplogroup exists
    """
    if os.path.isfile(input_text):
        try:
            text_result.insert(tk.END, "Raw DNA file detected! Processing...\n")
            root.update()

            user_snps = parse_user_dna(input_text)
            text_result.insert(tk.END, f"Extracted {len(user_snps)} valid mtDNA SNPs. Inferring haplogroup...\n")
            root.update()

            top_matches = infer_best_haplogroup(user_snps, tree_dict)

            text_result.insert(tk.END, "-" * 40 + "\n")
            text_result.insert(tk.END, "Haplogroup Inference Results:\n")
            for i, match in enumerate(top_matches):
                text_result.insert(tk.END, f" Top {i+1}: {match['haplogroup']} (Confidence: {match['score']:.1%})\n")
            text_result.insert(tk.END, "-" * 40 + "\n")

            if top_matches:
                user_hg = top_matches[0]['haplogroup']
                text_result.insert(tk.END, f"\nAuto-selected best match: 【{user_hg}】\n\n")
            else:
                text_result.insert(tk.END, "Could not infer haplogroup from the provided file.\n")
                return

        except Exception as e:
            messagebox.showerror("File read error", f"can't read file:\n{e}")
            return
            
    else:
        user_hg = input_text
        text_result.insert(tk.END, f"Manually entered haplogroup detected 【{user_hg}】\n\n")


    # Check if haplogroups exist in the tree
    if user_hg not in valid_hgs:
        text_result.insert(tk.END, f"No haplogroup found '{user_hg}'。\n")
        return

    
    """
    To find the closest ancient individuals:
    1. Iterate through the ancient database (aadr_df) row by row: Extract the target haplogroup for each ancient individual.
    2. Calculate the genetic distance using networkx (nx.shortest_path_length)
        Find the shortest path in the undirected graph (G) from user_hg (source) to ancient_hg (target)
        - up-forward, to tract the Most Recent Common Ancestor(in mtDNAtree17 graph)
        - then down-forward, to trace the ancient lineage(in aadr_df)
    3. Store the results:
        Save the calculated distances and ancient info into a list, convert to a DataFrame, sort by distance (ascending), and select the Top 5.
    4. Format the output string
    """
    text_result.insert(tk.END, f"Searching for ancient relatives closest to 【{user_hg}】...\n\n")
    root.update()

    distances = []
    for _, row in aadr_df.iterrows():
        ancient_hg = row['Haplogroup']
        try:
            dist = nx.shortest_path_length(G, source=user_hg, target=ancient_hg)
            distances.append({
                'Sample_ID': row['Sample_ID'],
                'Haplogroup': ancient_hg,
                'Distance': dist,
                'Location': row['Location'],
                'Age_BP': row['Age_BP']
            })
        except nx.NetworkXNoPath:  # if can't calculate distance, skip it
            continue
            
    results_df = pd.DataFrame(distances).sort_values(by='Distance').head(5)
    print(f"DEBUG: Found {len(results_df)} results")
    print(results_df[['Sample_ID', 'Distance']])
    
    text_result.insert(tk.END, "Match complete! Your Top 5 ancient relatives are as follows:\n")
    text_result.insert(tk.END, "-" * 50 + "\n")
    
    for i, (_, match) in enumerate(results_df.iterrows(), 1):
        result_str = (
            f"Top {i}: Sample id: {match['Sample_ID']} | Genetic distance: {match['Distance']} steps\n"
            f"      Haplogroup: {match['Haplogroup']}\n"
            f"      Period of life: Approximately {match['Age_BP']} years ago\n"
            f"      Location of origin: {match['Location']}\n"
        )+ "-" * 50 + "\n"
        print(f"DEBUG: {result_str}")
        text_result.insert(tk.END, result_str)


def main():
    global root, entry_hg, text_result

    # Draw GUI interface
    root = tk.Tk()
    root.title("mtDNA haplogroup matches system")
    root.geometry("600x650")

    label_inst = tk.Label(root, text="Please enter the haplogroup (e.g., X2c2) \nor the path to the results file. (e.g., C:\\user.txt):", font=("Arial", 11))
    label_inst.pack(pady=10)

    global entry_hg, text_result # clare global parameters
    entry_hg = tk.Entry(root, font=("Arial", 12), width=45)
    entry_hg.pack(pady=5)

    btn_search = tk.Button(root, text="Start matching...", font=("Arial", 12, "bold"), bg="lightblue", command=run_match)
    btn_search.pack(pady=10)

    text_result = tk.Text(root, font=("Courier", 10), width=65, height=16)
    text_result.pack(pady=10)

    root.mainloop()


if __name__ == "__main__":
    main()