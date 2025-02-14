import argparse
import os
import multiprocessing
import pandas as pd
from bin.utils import genbank_to_protein_fasta, make_diamond_db, run_diamond, make_presence_absence_table, calculate_jaccard_similarity_scores, construct_matrix, find_medoid, calculate_representation_scores

description = """
Find the most representative genome from a set of similar genomes, based on the presence/absence of genes.

Developer: Jordi Weijers
Affiliation: Institute of Biology Leiden, Leiden University
Please contact Jordi at s2228483@vuw.leidenuniv.nl if you have any issues.
Version: v0.2.0"""

def parse_arguments():
	""" Setup argument parsing and return parsed arguments."""
	parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
	input_group = parser.add_argument_group("Input Options")
	input_group.add_argument('-i', '--input_dir', type=str, required=True, help="Path to the input directory containing GenBank files.", metavar="")
	output_group = parser.add_argument_group("Output Options")
	output_group.add_argument('-o', '--output_dir', type=str, required=True, help="Path to save the output directory to.", metavar="")
	parser.add_argument('-qc', '--query_cover', type=float, default=70.0, help="Minimum percent query cover for diamond search. (Default: 70)", metavar="")
	parser.add_argument('-id', '--identity', type=float, default=35.0, help="Minimum percent identity for diamond search. (Default: 35)", metavar="")
	return parser.parse_args()

def validate_arguments(genbank_dir: str, output_dir: str, query_cover: float, identity: float) -> None:
	""" Validate the arguments."""
	if not os.path.isdir(genbank_dir):
	        raise ValueError(f"Error: The input directory '{genbank_dir}' does not exists. Please provide a directory containing GenBank files.")
	for filename in os.listdir(genbank_dir):
        	if not filename.endswith(".gbff"):
                	raise ValueError(f"Error: The input directory '{genbank_dir}' contains non-GenBank files. Please provide a directory ONLY containing GenBank files.")
	if os.path.exists(output_dir):
        	raise ValueError(f"Error: The output directory '{output_dir}' already exists. Please provide a new directory.")
	else:
        	os.makedirs(output_dir)
	if not (0 <= query_cover <= 100):
		raise ValueError("Error: --query_cover must be between 0 and 100.")

def process_genbank(filename, genbank_dir, pa_dir, blast_output_tsv, pa_tables):
	"""
	Process a GenBank file and make a presence absence table and calculate the representation score
	"""
	# Making the presence absence table
	basename, ext = os.path.splitext(filename)
	pa_table_path = os.path.join(pa_dir, f"{basename}.csv")
	pa_table = make_presence_absence_table(filename, blast_output_tsv)
	pa_table.to_csv(pa_table_path)
	pa_tables[basename] = pa_table

def main ():
	# Parse and Validate arguments
	args = parse_arguments()
	genbank_dir = args.input_dir
	output_dir = args.output_dir
	query_cover = args.query_cover
	identity = args.identity
	validate_arguments(genbank_dir, output_dir, query_cover, identity)

	# Initialize variables
	pa_tables = {}

	# Make protein FASTA files
	print("Making protein FASTA files")
	for filename in os.listdir(genbank_dir):
		genbank_file = os.path.join(genbank_dir, filename)
		basename, ext = os.path.splitext(filename)
		protein_fasta_dir = os.path.join(output_dir, "protein_fasta")
		os.makedirs(protein_fasta_dir, exist_ok=True)
		protein_fasta = os.path.join(protein_fasta_dir, f"{basename}.fna")
		genbank_to_protein_fasta(genbank_file, protein_fasta)

	# Make a combined protein FASTA file
	all_protein_fasta = os.path.join(output_dir, "all_protein.fna")
	with open(all_protein_fasta, "w") as outfile:
		for fasta_file in os.listdir(protein_fasta_dir):
			if fasta_file.endswith(".fna"):
				fasta_path = os.path.join(protein_fasta_dir, fasta_file)
				with open(fasta_path, "r") as infile:
					outfile.write(infile.read())
	# Make a DIAMOND database
	print("Making the DIAMOND database")
	diamond_db = os.path.join(output_dir, "all_protein.dmnd")
	make_diamond_db(all_protein_fasta, diamond_db)

	# Run DIAMOND all vs all
	print("Running the DIAMOND all vs all")
	blast_output_tsv = os.path.join(output_dir, "blast.tsv")
	run_diamond(all_protein_fasta, diamond_db, blast_output_tsv, query_cover, identity)

	# Make a presence absence table for each GenBank file based on the DIAMOND output
	print("Making the presence-absence tables")
	pa_dir = os.path.join(output_dir, "presence_absence")
	os.makedirs(pa_dir, exist_ok=True)
	# Use multiprocessing for parallel processing
	manager = multiprocessing.Manager()
	pa_tables = manager.dict()
	pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
	tasks = [(filename, genbank_dir, pa_dir, blast_output_tsv, pa_tables) for filename in os.listdir(genbank_dir)]
	pool.starmap(process_genbank, tasks)
	pool.close()
	pool.join()
	pa_tables = dict(pa_tables)

	# Calculating the Jaccard similarity scores
	print("Calculating the Jaccard similarity scores")
	jaccard_similarity_scores = calculate_jaccard_similarity_scores(pa_tables)
	# Write the Jaccard similarity scores to matrix files
	jaccard_matrix_dir = os.path.join(output_dir, "jaccard_similarity_matrices")
	os.makedirs(jaccard_matrix_dir, exist_ok=True)
	for ref_genbank_filename, jac_sim_scores in jaccard_similarity_scores.items():
		jac_sim_df = pd.DataFrame(list(jac_sim_scores.items()), columns=["Other_GenBank", "Jaccad_Similarity_Score"])
		jac_sim_df.insert(0, "Reference_GenBank", ref_genbank_filename)
		output_file = os.path.join(jaccard_matrix_dir, f"{ref_genbank_filename}_jaccard_matrix.csv")
		jac_sim_df.to_csv(output_file, index=False)
	# Concatenate the Jaccard matrices
	jaccard_matrix = construct_matrix(jaccard_similarity_scores)
	# Write the concatenated Jaccard matrix to file
	jaccard_matrix_path = os.path.join(output_dir, "jaccard_matrix.csv")
	jaccard_matrix.to_csv(jaccard_matrix_path)
	# Calculate the representation scores
	representation_scores_df = calculate_representation_scores(jaccard_matrix)
	output_file = os.path.join(output_dir, "representation_scores.tsv")
	representation_scores_df.to_csv(output_file, sep="\t", index=False)
	representative = find_medoid(jaccard_matrix)

	print(f"All steps completed. Results saved to '{output_dir}'.")
	print(f"The representative is: {representative}")

if __name__ == "__main__":
	main()
