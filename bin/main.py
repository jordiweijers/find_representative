import argparse
import os
from find_representative.utils import genbank_to_protein_fasta, make_diamond_db, run_diamond, make_presence_absence_table, calculate_representation_score

def parse_arguments():
	""" Setup argument parsing and return parsed arguments."""
	parser = argparse.ArgumentParser(description="Extract protein sequences from a GenBank file and save them to a FASTA file.")
	parser.add_argument('--input_dir', type=str, required=True, help="Path to the input directory containing GenBank files.")
	parser.add_argument('--output_dir', type=str, required=True, help="Path to save the output directory to.")
	args = parser.parse_args()
	return parser.parse_args()

def validate_arguments(genbank_dir: str, output_dir: str) -> None:
	""" Validate the arguments."""
	if not os.path.isdir(genbank_dir):
	        raise ValueError(f"The input directory '{genbank_dir}' does not exists. Please provide a directory containing GenBank files.")
	for filename in os.listdir(genbank_dir):
        	if not filename.endswith(".gbff"):
                	raise ValueError(f"The input directory '{genbank_dir}' contains non-GenBank files. Please provide a directory ONLY containing GenBank files.")
	if os.path.exists(output_dir):
        	raise ValueError(f"The output directory '{output_dir}' already exists. Please provide a new directory.")
	else:
        	os.makedirs(output_dir)

def main ():
	# Parse and Validate arguments
	args = parse_arguments()
	genbank_dir = args.input_dir
	output_dir = args.output_dir
	validate_arguments(genbank_dir, output_dir)

	# Initialize variables
	rep_score_dict = {}

	# Make protein FASTA files
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
	diamond_db = os.path.join(output_dir, "all_protein.dmnd")
	make_diamond_db(all_protein_fasta, diamond_db)

	# Run DIAMOND all vs all
	blast_output_tsv = os.path.join(output_dir, "blast.tsv")
	run_diamond(all_protein_fasta, diamond_db, blast_output_tsv)

	# Make a presence absence table for each GenBank file based on the DIAMOND output
	pa_dir = os.path.join(output_dir, "presence_absence")
	os.makedirs(pa_dir, exist_ok=True)
	for filename in os.listdir(genbank_dir):
		basename, ext = os.path.splitext(filename)
		pa_table_path = os.path.join(pa_dir, f"{basename}.csv")
		pa_table = make_presence_absence_table(filename, blast_output_tsv)
		pa_table.to_csv(pa_table_path)
		# Calculate the representation score
		rep_score = calculate_representation_score(filename, pa_table)
		rep_score_dict[basename] = rep_score
	# Sort the representation dictionary in descending order
	sorted_rep_score_dict = sorted(rep_score_dict.items(), key=lambda x: x[1], reverse=True)
	# Print sorted representation dictionary
	print("\nRepresentation Scores:")
	for genbank_filename, rep_score in sorted_rep_score_dict:
		print(f"{genbank_filename}: {rep_score:.2f}")

if __name__ == "__main__":
	main()
