#!/usr/bin/env python3

"""
Author: Yaqi Jiao
Date: 12th March, 2026

user_data_process.py
----------------
Description:
    This module handles the ingestion and preprocessing of raw, consumer-grade DNA 
    sequencing files (e.g., 23andMe, AncestryDNA formats). Its primary function is 
    to extract a clean, high-quality set of mitochondrial DNA (mtDNA) single nucleotide 
    polymorphisms (SNPs) for downstream haplogroup inference.

    Key functionalities:
    1. Data Parsing: Efficiently reads tab-separated raw genotype files using pandas, 
    automatically bypassing metadata headers (lines starting with '#').
    2. mtDNA Isolation: Standardizes and filters mitochondrial data across various 
    chromosome naming conventions used by different sequencing companies (e.g., 'MT', 'M', '26').
    3. Quality Control (QC): Cleans the genotype data by filtering out "no-calls" 
    (e.g., '--', '0', 'NaN') and invalid inputs, ensuring only standard nucleotides 
    (A, T, C, G) are processed.
    4. Data Structuring: Outputs an optimized Python dictionary mapping genomic positions 
    to their corresponding alleles, serving as the direct input for the inference engine.


"""

import pandas as pd

def parse_user_dna(file_path):
    """
    Reads the user's raw DNA sequencing file and extracts the mutation sites in the mtDNA.
    Returns a dictionary: {position: 'allele'}, for example, {263: 'G', 752566: 'A'}
    """
    print(f"Loading user DNA data from: {file_path} ...")
    
    column_names = ['rsid', 'chromosome', 'position', 'genotype']
    df = pd.read_csv(
        file_path, 
        sep='\t', 
        comment='#',          # ignore lines starts with "#"
        names=column_names,
        low_memory=False
    )

    df['chromosome'] = df['chromosome'].astype(str).str.upper()
    mt_df = df[df['chromosome'].isin(['MT', 'M', '26'])]  # extract mtDNA data
    user_mt_snps = {}

    for _, row in mt_df.iterrows():
        try: 
            # extract position and genotype information
            pos = int(row['position'])
            geno = str(row['genotype']).strip().upper()
            
            # quality control
            if len(geno) > 0 and geno[0] in ['A', 'T', 'C', 'G']: 
                allele = geno[0]
                user_mt_snps[pos] = allele
                
        except ValueError:
            # ignore abnormal values
            continue
    
    print(f"Success! Extracted {len(user_mt_snps)} valid mtDNA positions from user data.")
    return user_mt_snps
