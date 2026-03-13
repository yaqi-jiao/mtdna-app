# mtDNA Matcher

**mtDNA Matcher** is a Python-based desktop application designed to infer a user's mitochondrial DNA (mtDNA) haplogroup from raw genotype data and find their ancient human relatives. Developed as part of the **BINP29** course at Lund University.

This tool utilizes ancient DNA databases (like AADR, AGDP) and the human mtDNA phylogenetic tree (Build 17) to accurately predict haplogroups from raw SNPs, calculate genetic distances, and identify the closest ancient matches to a user's genetic profile.

## Features
* **User-Friendly GUI**: Built with Tkinter for an intuitive desktop experience.
* **Advanced Haplogroup Inference**: Automatically extracts mtDNA-specific SNPs from user-provided raw genotype files (e.g., 23andMe raw data) and infers the most likely haplogroup by mapping mutations against the Build 17 tree.
* **Flexible Input Modes**: 
  * Direct string input (e.g., `X2c2`).
  * Raw data parsing (supports `.txt` files containing unphased genotype data).
* **Phylogenetic Network Analysis**: Uses `NetworkX` to calculate true evolutionary distances based on the tree topology, rather than relying on simple string matching.
* **Automated Data Processing**: Includes standalone scripts to parse raw ancient metadata and construct complex tree topologies from scratch.

## Repository Structure
```text
mtDNA_Project/
├── pyproject.toml              # Modern Python package configuration
├── README.md                   # Project documentation
├── Test/                       # Test files for users
│   ├── Test1.txt
│   ├── Test2.txt         
│   └── Test3.txt               
└── mtDNAmatcher/
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
git clone [https://github.com/yaqi-jiao/mtdna-app.git](https://github.com/yaqi-jiao/mtdna-app.git) 
cd mtDNA_Project
```

Create and activate a virtual environment (optional but recommended):

```bash
python -m venv test
# On Windows:
test\Scripts\activate
# On macOS/Linux:
source test/bin/activate
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

