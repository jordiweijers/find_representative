##Find Representative

Find Representative is a tool used to select the most representative genome from a set of similar genomes. It does this based on the presence/absence of genes and gives a higher score when the genes in the genome are present across the whole set of genomes.

##Installation
Installation can be performed via cond

```bash
# 1. clone the repository
git clone https://github.com/jordiweijers/find_representative.git
cd find_representative

# 2. create conda environment using yaml file and activate it.
conda env create -f find_representative_env.yml
conda activate find_representative

# 3. install the python package
pip install .
```

##Usage

Find Representative takes a directory containing genbank files and outputs a directory with the output files

```
optional arguments:
  -h, --help            show this help message and exit
  -qc , --query_cover   Minimum percent query cover for diamond search. (Default: 70)
  -id , --identity      Minimum percent identity for diamond search. (Default: 35)

Input Options:
  -i , --input_dir      Path to the input directory containing GenBank files.

Output Options:
  -o , --output_dir     Path to save the output directory to.
```
