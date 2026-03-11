#!/usr/bin/env python3

"""
Author: Yaqi Jiao
Date: 9th March, 2026

mtDNAmatcher.py
----------------

Description: This application identifies the closest ancient human relatives based on a user's mtDNA haplogroup.
    1. Load background data: Reconstructs the mtDNA phylogenetic tree as an undirected graph and loads the ancient DNA database (AADR).
    2. Read user input: Acquires input from the GUI text entry.
        - If user input is a file path: Opens the file and extracts the haplogroup string following the "mtDNA: " pattern using regex.
        - If user input is a standard string: Treats the input directly as the target haplogroup.
    3. Calculate and Match: Validates the haplogroup, calculates the genetic distance (shortest path in the tree graph) between the user and ancient individuals, and sorts them to find the closest matches.

Input:
    Via GUI text entry. Accepts either:
    - A string representing an mtDNA haplogroup (e.g., "X2c2").
    - A valid file path to a text file containing the user's haplogroup in the format "mtDNA: [haplogroup]".

Output:
    Displays the Top 5 closest ancient matches in the GUI text area, including their Sample ID, Genetic distance (steps in the phylogenetic tree), Haplogroup, Period of life (Age BP), and Location of origin.

Usage Example:
    $ python mtDNAmatcher.py
    # A GUI window will launch. 
    # Enter "X2" or "C:\\user_results.txt" into the input box and click "Start matching...".
"""

import pandas as pd
import networkx as nx
import tkinter as tk
from tkinter import messagebox
import os
import re


# Load background data
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
print("Loading data...")
aadr_file = os.path.join(BASE_DIR, "Data", "Clean_metadata.tsv") 
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
        if the file exists, use re to extract content after "mtDNA:" as hg_name
    2. If the input is not a file path, then use input as hg_name
    3. Search valid_hgs with user_hg, to make sure the haplogroup exists
    """
    if os.path.isfile(input_text):
        try:
            with open(input_text, 'r', encoding='utf-8') as f:
                content = f.read()
                match = re.search(r'mtDNA:\s*(\S+)', content, re.IGNORECASE)
                if match:
                    user_hg = match.group(1)
                    text_result.insert(tk.END, f"file input identified, extract mtDNA haplogroup: 【{user_hg}】\n\n")
                else:
                    text_result.insert(tk.END, "Can find 'mtDNA: XXX' format conten, check file format\n")
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

    text_result.insert(tk.END, f"Looking for the ancient people closest to the haplogroup 【{user_hg}】...\n\n")
    
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
    # Draw GUI interface
    root = tk.Tk()
    root.title("mtDNA haplogroup matches system")
    root.geometry("550x500")

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