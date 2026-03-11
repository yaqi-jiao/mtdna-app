# mtDNA Matcher

**mtDNA Matcher** is a Python-based desktop application designed to find ancient human relatives based on mitochondrial DNA (mtDNA) haplogroups. Developed as part of the **BINP29** course at Lund University.

This tool utilizes ancient DNA databases (like AGDP) and the human mtDNA phylogenetic tree (Build 17) to calculate genetic distances and identify the closest ancient matches to a user-provided haplogroup.

## Features
* **User-Friendly GUI**: Built with Tkinter for an intuitive desktop experience.
* **Dual Input Modes**: 
  * Direct string input (e.g., `X2c2`).
  * File parsing (supports `.txt` files containing haplogroup information).
* **Phylogenetic Network Analysis**: Uses `NetworkX` to traverse the mtDNA tree and calculate true evolutionary distances, rather than just string matching.
* **Automated Data Cleaning**: Includes standalone scripts to parse and clean raw metadata and tree topologies.

## Repository Structure
```text
mtDNA_Project/
├── pyproject.toml              # Modern Python package configuration
├── README.md                   # Project documentation
├── example_data/               # Test files for users
│   └── test_sample.txt         # Example input file (contains: mtDNA: X2)
└── mtDNAmatcher/               # Main source code and databases
    ├── __init__.py
    ├── mtDNAmatcher.py         # Core GUI and matching algorithm
    ├── metadata_parser.py      # Script to clean raw AGDP metadata
    ├── mtDNA_tree17_build.py   # Script to parse the mtDNA Build 17 tree
    └── Data/                   # Cleaned TSV datasets required for the app
```
## Getting started

### 1. Prerequisites

  - python 3.8+
  - It is highly recommended to use a python environment

### 2. Installation

Clone this repository to your local machine:

```bash
git clone 
cd mtDNA_Project
```

Install the application and its dependencies (pandas, networkx, openpyxl) via pip:

```bash
# Make sure your virtual environment is activated before running this command
pip install .
```

## Usage

Once installed, the application is registered to your system path. You can launch it from any directory in your terminal by simply typing:

```bash
mtdna-app
```

### Testing the Application

When the graphical interface opens, you can test it in two ways:

1. **Direct Input**: Type a haplogroup like X2 or H1a into the text box and click "Start matching...".
2. **File Input:** Type the relative path to the provided test file to simulate reading from a sequencing result:
   ```bash
  TestData/test1.txt
   ```

## Data processing

If you want to update the background databases, you can run the parsing scripts directly:

```bash
python mtDNAmatcher/metadata_parser.py
python mtDNAmatcher/mtDNA_tree17_build.py
```

