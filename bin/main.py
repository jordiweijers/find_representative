import argparse
import os
import multiprocessing
from bin.utils import genbank_to_protein_fasta, make_diamond_db, run_diamond, make_presence_absence_table, calculate_representation_score

def parse_arguments():
	""" Setup argument parsing and return parsed arguments."""
	parser = argparse.ArgumentParser(description="Find the most representative genome from a set of similar genomes, based on the presence/absence of genes")
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

def process_genbank(filename, genbank_dir, pa_dir, blast_output_tsv, rep_score_dict):
	"""
	Process a GenBank file and make a presence absence table and calculate the representation score
	"""
	# Making the presence absence table
	basename, ext = os.path.splitext(filename)
	pa_table_path = os.path.join(pa_dir, f"{basename}.csv")
	pa_table = make_presence_absence_table(filename, blast_output_tsv)
	pa_table.to_csv(pa_table_path)
	# Calculate the representation score
	rep_score = calculate_representation_score(filename, pa_table)
	rep_score_dict[basename] = rep_score

def main ():
	# Parse and Validate arguments
	args = parse_arguments()
	genbank_dir = args.input_dir
	output_dir = args.output_dir
	query_cover = args.query_cover
	identity = args.identity
	validate_arguments(genbank_dir, output_dir, query_cover, identity)

	# Initialize variables
	rep_score_dict = {}

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
	rep_score_dict = manager.dict()
	pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
	tasks = [(filename, genbank_dir, pa_dir, blast_output_tsv, rep_score_dict) for filename in os.listdir(genbank_dir)]
	pool.starmap(process_genbank, tasks)
	pool.close()
	pool.join()
	rep_score_dict = dict(rep_score_dict)
	# Sorting the representation dictionary from highest to lowest
	sorted_rep_score_dict = sorted(rep_score_dict.items(), key=lambda x: x[1], reverse=True)
	# Print and wirte sorted representation dictionary to a file
	output_filename = os.path.join(output_dir, "representation_scores.txt")
	with open(output_filename, "w") as file:
		print("\nRepresentation Scores:")
		file.write("Representation Scores:\n")
		for genbank_filename, rep_score in sorted_rep_score_dict:
			output_line = f"{genbank_filename}: {rep_score:.2f}"
			print(output_line)
			file.write(output_line + "\n")
	print(f"Results saved to {output_filename}")

if __name__ == "__main__":
	main()
