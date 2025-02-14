import os
from Bio import SeqIO
import shutil
import subprocess
import pandas as pd
import numpy as np
from sklearn.manifold import MDS

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
	basename =  os.path.splitext(os.path.basename(genbank_file))[0]
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
	basename = os.path.splitext(genbank_filename)[0]
	columns = ["qseqid", "sseqid"]
	blast_df = pd.read_csv(blast_output, sep="\t", header=None, names=columns)
	blast_df = blast_df[blast_df["qseqid"].str.endswith(basename)]
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

def calculate_jaccard_similarity_scores(pa_tables: dict[str, pd.DataFrame]) -> dict[str, dict[str, float]]:
	"""
	Calculate the Jaccard similarity scores based on a dictionary containing all the presence-absence tables.
	Args:
		pa_tables: A Dictionary containing all the presence-absence tables.
	Returns:
		jaccard_similarity_scores: A Dictionary containing the different Jaccard similarity matrices for each presence-absence table.
	"""
	jaccard_similarity_scores = {}
	for ref_genbank_filename, pa_table in pa_tables.items():
		if ref_genbank_filename not in pa_table.index:
			raise ValueError(f"Reference genome '{ref_genbank_filename}' not found in the pa_table index '{pa_table}'.")
		gene_frequency = {}
		for gene_name in pa_table.columns:
			gene_freq = (pa_table[gene_name] > 0).sum()
			gene_frequency[gene_name] = gene_freq
		max_frequency = max(gene_frequency.values())
		weights = {}
		for gene_name, freq in gene_frequency.items():
			weights[gene_name] = 1 - (freq / max_frequency)
		ref_vector = pa_table.loc[ref_genbank_filename]
		jac_sim_scores = {}
		for other_genbank_filename in pa_table.index:
			other_vector = pa_table.loc[other_genbank_filename]
			weighted_intersection = 0
			weighted_union = 0
			for gene_name in pa_table.columns:
				weight = weights[gene_name]
				count_ref = ref_vector[gene_name]
				count_other = other_vector[gene_name]
				weighted_intersection += weight * min(count_ref, count_other)
				weighted_union += weight * max(count_ref, count_other)
			jac_sim = weighted_intersection / weighted_union if weighted_union != 0 else 0.0
			jac_sim_scores[other_genbank_filename] = jac_sim
		jaccard_similarity_scores[ref_genbank_filename] = jac_sim_scores
	return jaccard_similarity_scores

def construct_matrix(jaccard_similarity_scores: dict[str, dict[str, float]]) -> pd.DataFrame:
	"""
	Construct a complete and symmetrical Jaccard similarity matrix by concatenating and averaging the the Jaccard scores.
	Args:
		jaccard_similarity_scores: A Dictionary containing the different Jaccard similarity matrices for each presence-absence table.
	Returns:
		jaccard_matrix: A Pandas DataFrame containing the complete Jaccard similarity matrix.
	"""
	genomes = list(jaccard_similarity_scores.keys())
	jaccard_matrix = pd.DataFrame(index=genomes, columns=genomes, dtype=float)
	for genome_i in genomes:
		for genome_j in genomes:
			if genome_i == genome_j:
				jaccard_matrix.loc[genome_i, genome_j] = jaccard_similarity_scores[genome_i][genome_j]
			else:
				score_i_j = jaccard_similarity_scores[genome_i].get(genome_j)
				score_j_i = jaccard_similarity_scores[genome_j].get(genome_i)
				if score_i_j is not None and score_j_i is not None:
					average_score = (score_i_j + score_j_i) / 2
				elif score_i_j is None:
					average_score = score_j_i
				elif score_j_i is None:
					average_score = score_i_j
				else:
					average_score = 0
				jaccard_matrix.loc[genome_i, genome_j] = average_score
	return jaccard_matrix

def find_medoid(jaccard_matrix: pd.DataFrame) -> str:
	"""
	Find the medoid directly from the jaccard_matrix, which is the point for which the total sum of distances is minimal, and therefore the most representative.
	Args:
		jaccard_matrix: A Pandas DataFrame containing the complete Jaccard similarity matrix.
	Returns:
		medoid_index: The index label of the medoid.
	"""
	distance_matrix = 1 - jaccard_matrix
	medoid_index = distance_matrix.sum(axis=1).idxmin()
	return medoid_index

def calculate_representation_scores(jaccard_matrix: pd.DataFrame) -> dict[str, float]:
	"""
	Calculate the distances from 1 point to all others which will be the representation scores. The lower the score the better.
	Args:
		jaccard_matrix: A Pandas DataFrame containing the complete Jaccard similarity matrix.
	Returns:
		representation_scores: A Pandas DataFrame containing the representation scores.
	"""
	distance_matrix = 1 - jaccard_matrix
	rep_scores = distance_matrix.sum(axis=1)
	representation_scores = pd.DataFrame({"Genome": rep_scores.index, "Representation_Score": rep_scores.values})
	representation_scores = representation_scores.sort_values(by="Representation_Score", ascending=True).reset_index(drop=True)
	return representation_scores
