import os
from Bio import SeqIO
import shutil
import subprocess
import pandas as pd

def genbank_to_protein_fasta(genbank_file: str, protein_fasta: str) -> None:
	"""
	Converts a GenBank file to a FASTA file containing protein coding sequences.

	Args:
		genbank_file: Path to input GenBank file.
		protein_fasta: Path to output protein FASTA file.
	Returns:
		None
	"""
	lines = []
	basename =  os.path.basename(genbank_file)
	for record in SeqIO.parse(genbank_file, "genbank"):
		for feature in record.features:
			if feature.type == "CDS":
				locus_tag = feature.qualifiers.get("locus_tag",  ["N/A"])[0].split("|")[0]
				location = str(feature.location)
				location = location.replace(" ", "")
				lines.append(">" + locus_tag + "|" + location + "|" + basename)
				protein_seq = str(feature.qualifiers.get("translation")[0])
				lines.append(protein_seq)
	with open(protein_fasta, "a") as file:
		for line in lines:
			file.write(line + '\n')
		file.close()

def make_diamond_db(protein_fasta: str, dmnd_db: str, dmnd_path: str = "diamond") -> None:
	"""
	Makes a DIAMOND databases from a FASTA file containing protein sequences.

	Args:
		protein_fasta: Path to input protein FASTA file.
		dmnd_db: Path to output DIAMOND database.
		dmnd_path: Path to the DIAMOND executable (default is 'diamond' if it is in the PATH).
	Returns:
		None
	"""

	if not shutil.which(dmnd_path):
		raise FileNotFoundError(f"DIAMOND executable not found at {dmnd_path}. Please ensure it is installed and in your PATH.")
	command = [
		dmnd_path, 'makedb',
		"--in", protein_fasta,
		"--db", dmnd_db,
		"--quiet"
	]
	subprocess.run(command, check=True)

def run_diamond(query: str, dmnd_db: str, tsv_file: str, query_cover: float, identity: float, dmnd_path: str = "diamond") -> None:
	"""
	Runs DAIMOND BLAST to align protein sequences against a DIAMOND database.

	Args:
		query: Path to a input query protein FASTA file.
		dmnd_db: Path to a input DIAMOND database.
		tsv_file: Path to tsv file containing the output.
		query_cover: A float representing the minimum query cover percentage.
		identity: A float representing the minimum identity percentage.
		dmnd_path: Path to the DIAMOND executable (default is 'diamond'if it is in the PATH).
	Returns:
		None
	"""
	if not shutil.which(dmnd_path):
		raise FileNotFoundError(f"DIAMOND executable not found at {dmnd_path}. Please ensure it is installed and in your PATH.")

	command = [
		dmnd_path, "blastp",
		"--query", query,
		"--db", dmnd_db,
		"--out", tsv_file,
		"--query-cover", str(query_cover),
		"--id", str(identity),
		"-k0",
		"--quiet",
		"--outfmt", "6", "qseqid", "sseqid"
	]
	subprocess.run(command, check=True)

def make_presence_absence_table(genbank_filename: str, blast_output: str) -> pd.DataFrame:
	"""
	Generate a presence absence table for a given GenBank filename based on the DIAMOND blast output.
	Args:
		genbank_filename: The name of a GenBank file.
		blast_output: Path to the DIAMOND BLAST output file.
	Returns:
		presence_absence: A Pandas DataFrame representing the presence absence table.
	"""
	columns = ["qseqid", "sseqid"]
	blast_df = pd.read_csv(blast_output, sep="\t", header=None, names=columns)
	blast_df = blast_df[blast_df["qseqid"].str.endswith(genbank_filename)]
	current_genes = blast_df["qseqid"].unique()
	other_genbank_names = blast_df["sseqid"].apply(lambda x: x.split("|")[-1]).unique()
	presence_absence = pd.DataFrame(0, index=other_genbank_names, columns=current_genes)
	for name in other_genbank_names:
		for qseqid in current_genes:
			hits = blast_df[
				(blast_df["qseqid"] == qseqid) &
				(blast_df["sseqid"].str.split("|").str[-1] == name)
			]
			presence_absence.at[name, qseqid] = len(hits)
	return presence_absence

def calculate_gene_family_base_value(gene_name: str, pa_table: pd.DataFrame) -> float:
	"""
	Calculate the gene family base for a given gene and presence-absence table.
	Args:
		gene_name: The header of a gene (column name in the pa_table).
		pa_table: A Pandas DataFrame containing the presence absence table.
	Returns:
		gene_family_base: A float representing the gene family base value (basically the conservation percentage).
	"""
	if gene_name not in pa_table.columns:
		raise ValueError(f"Gene '{gene_name}' not found in the presence-absence table")
	presence_count = (pa_table[gene_name] > 0).sum()
	total_count = len(pa_table)
	base = presence_count / total_count
	return base

def calculate_gene_family_modifier(gene_name: str, genbank_filename: str, pa_table: pd.DataFrame) -> float:
	"""
	Calculate the gene family modifier for a specific genbank filename based on the deviation from the most common (mode) copy number in the presence-absence table.
        Args:
                gene_name: The header of a gene (column name in the pa_table).
		genbank_filename: The name for which GenBank file the gene family modifier needs to be calculated (Normally this will be the same as the inputted GenBank filename for the presence-absence table).
                pa_table: A Pandas DataFrame containing the presence absence table.
	Returns:
		gene_family_modifier: A float representing the gene family modifier value.
	"""
	if gene_name not in pa_table.columns:
		raise ValueError(f"Gene '{gene_name}' not found in the presence-absence table")
	filtered_values = pa_table[gene_name][pa_table[gene_name] != 0]
	mode = filtered_values.mode().iloc[0]
	frequency = pa_table.loc[genbank_filename, gene_name]
	modifier = 1 - (abs(frequency - mode)/((frequency + mode) / 2))
	return modifier

def calculate_representation_score(genbank_filename: str, pa_table: str) -> float:
	"""
	Calculate the final score to make the best hit.
	Args:
		genbank_filename: The name for which GenBank file the representation score needs to be calculated.
		pa_table: A Pandas DataFrame containing the presence absence table.
	Returns:
		score: A float representing a score that says something about how representative a given GenBank file is compared to others.
	"""
	if genbank_filename not in pa_table.index:
		raise ValueError(f"GenBank file '{genbank_filename}' not found in the presence-absence table")
	score = 0
	for gene_name in pa_table.columns:
		base = calculate_gene_family_base_value(gene_name, pa_table)
		modifier = calculate_gene_family_modifier(gene_name, genbank_filename, pa_table)
		score += base * modifier
	return score
