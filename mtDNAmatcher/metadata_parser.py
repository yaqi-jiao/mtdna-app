#!/usr/bin/env python3

"""
Author: Yaqi Jiao
Date: 9th March, 2026

metadata_parser.py
----------------

Description: 
    This script parses and cleans the raw ancient DNA metadata file (AGDP dataset). 
    It performs the following operations:
    1. Loads the raw Excel file.
    2. Extracts essential columns required for the matching system: sample ID, age, locality, and mtDNA haplogroup.
    3. Renames the lengthy original column headers into standardized, easy-to-use variables (Sample_ID, Age_BP, Location, Haplogroup).
    4. Filters out any records that lack a valid mtDNA haplogroup (removing NaNs, "n/a", "..", "-", etc.).
    5. Exports the cleaned dataset as a TSV file for downstream phylogenetic analysis.

Input:
    Raw metadata in Excel format (e.g., AGDP.metadata.xlsx).

Output:
    A cleaned, tab-separated values file (Clean_metadata.tsv) containing only samples with valid mtDNA haplogroups.

Usage Example:
    $ python metadata_parser.py
    # The script will output the number of extracted samples and the path to the saved TSV file.

"""

import pandas as pd
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.join(BASE_DIR, "Data", "AGDP.metadata.xlsx")
output_file = os.path.join(BASE_DIR, "Data", "Clean_metadata.tsv")

print(f"Loading data: {input_file} ...")
df = pd.read_excel(input_file)

target_columns = [
    'I-ID', 
    'Date mean in BP [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]', 
    'Locality', 
    'mtDNA haplogroup'
]

df = df[target_columns].rename(columns={
    'I-ID': 'Sample_ID',
    'Date mean in BP [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]': 'Age_BP',
    'Locality': 'Location',
    'mtDNA haplogroup': 'Haplogroup'
})

df = df.dropna(subset=['Haplogroup'])
df = df[~df['Haplogroup'].isin(['..', 'n/a', 'N/A', '-', ''])]

df.to_csv(output_file, sep='\t', index=False)

print(f"Process finished!  {len(df)} samples are extracted")
print(f"File saved as: {output_file}")
