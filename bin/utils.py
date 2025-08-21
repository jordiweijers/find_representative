import logging
import os
from Bio import SeqIO
import shutil
import subprocess
import csv
import pandas as pd
import numpy as np
from sklearn.manifold import MDS

logger = logging.getLogger(__name__)

def genbank_to_protein_fasta(genbank_file: str, protein_fasta: str, log_file: str) -> None:
	"""
	Converts a GenBank file to a FASTA file containing protein coding sequences.
	Skips CDS features without translations (e.g., pseudogenes or partial sequences).
	Logs summary statistics to a specified log file.

	Args:
		genbank_file: Path to input GenBank file.
		protein_fasta: Path to output protein FASTA file.
		log_file: Path to log file to append processing statistics.
	"""
	if not os.path.exists(genbank_file):
		logger.error(f"GenBank file not found: {genbank_file}.")
		raise FileNotFoundError(f"GenBank file not found: {genbank_file}.")

	lines = []
	basename = os.path.splitext(os.path.basename(genbank_file))[0]
	written_count = 0
	skipped_no_translation = 0
	skipped_pseudogene = 0

	try:
		for record in SeqIO.parse(genbank_file, "genbank"):
			for feature in record.features:
				if feature.type != "CDS":
					continue

				# Skip pseudogenes
				if "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers:
					skipped_pseudogene += 1
					continue

				# Get translation
				translation = feature.qualifiers.get("translation")
				if not translation:
					skipped_no_translation += 1
					continue

				locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0].split("|")[0]
				location = str(feature.location).replace(" ", "")
				protein_seq = str(translation[0]).replace(" ", "").replace("\n", "")

				lines.append(f">{locus_tag}|{location}|{basename}")
				lines.append(protein_seq)
				written_count += 1

		# Append sequences to output FASTA
		with open(protein_fasta, "a") as file:
			file.write("\n".join(lines) + "\n")

		# Append log entry
		with open(log_file, "a") as log:
			log.write(
				f"{basename}: Wrote {written_count} protein sequences | "
				f"Skipped {skipped_no_translation} with no translation | "
				f"Skipped {skipped_pseudogene} pseudogenes\n"
			)
	except Exception as e:
		logger.error(f"Error while processing {genbank_file}: {e}.", exc_info=True)
		raise

def make_diamond_db(protein_fasta: str, dmnd_db: str, cpus: int, dmnd_path: str = "diamond") -> None:
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
		logger.error(f"DIAMOND executable not found at {dmnd_path}. Please ensure it is installed an in your PATH.")
		raise FileNotFoundError(f"DIAMOND executable not found at {dmnd_path}. Please ensure it is installed and in your PATH.")
	command = [
		dmnd_path, 'makedb',
		"--in", protein_fasta,
		"--db", dmnd_db,
		"--threads", str(cpus),
		"--quiet"
	]
	try:
		subprocess.run(command, check=True)
	except subprocess.CalledProcessError as e:
		logger.error(f"DIAMOND BLAST failed: {e}.")
		raise RuntimeError(f"DIAMOND BLAST failed: {e}.")

def run_diamond(query: str, dmnd_db: str, tsv_file: str, query_cover: float, identity: float, cpus: int, dmnd_path: str = "diamond") -> None:
	"""
	Runs DAIMOND BLAST to align protein sequences against a DIAMOND database.

	Args:
		query: Path to a input query protein FASTA file.
		dmnd_db: Path to a input DIAMOND database.
		tsv_file: Path to tsv file containing the output.
		query_cover: A float representing the minimum query cover percentage.
		identity: A float representing the minimum identity percentage.
		cpues: The Number of CPUs to use.
		dmnd_path: Path to the DIAMOND executable (default is 'diamond'if it is in the PATH).
	Returns:
		None
	"""
	if not shutil.which(dmnd_path):
		logger.error(f"DIAMOND executable not found at {dmnd_path}. Please ensure it is installed an in your PATH.")
		raise FileNotFoundError(f"DIAMOND executable not found at {dmnd_path}. Please ensure it is installed and in your PATH.")

	command = [
		dmnd_path, "blastp",
		"--query", query,
		"--db", dmnd_db,
		"--out", tsv_file,
		"--query-cover", str(query_cover),
		"--id", str(identity),
		"--threads", str(cpus),
		"-k0",
		"--quiet",
		"--outfmt", "6", "qseqid", "sseqid"
	]
	try:
		subprocess.run(command, check=True)
	except subprocess.CalledProcessError as e:
		logger.error(f"DIAMOND BLAST failed: {e}.")
		raise RuntimeError(f"DIAMOND BLAST failed: {e}.")

def make_presence_absence_table(genbank_filename: str, blast_output: str) -> tuple[np.ndarray, dict, dict]:
	"""
	Generate a presence absence table for a given GenBank filename based on the DIAMOND blast output.
	Args:
		genbank_filename: The name of a GenBank file.
		blast_output: Path to the DIAMOND BLAST output file.
	Returns:
		presence_absence: NumPy array of shape (num_s_genomes, num_q_genes)
		genome_map_index: A Dictionary mapping the genome name to the row index
		gene_map_index: A Dictionary mapping the gene name to the column index
	"""
	if not os.path.exists(blast_output):
		logger.error(f"BLAST TSV file not found: {blast_output}.")
		raise FileNotFoundError(f"BLAST TSV file not found: {blast_output}.")

	basename = os.path.splitext(genbank_filename)[0]
	genes = []
	genomes = []
	with open(blast_output, newline = '',) as file:
		reader = csv.reader(file, delimiter = "\t")
		for row in reader:
			qseqid, sseqid = row[0], row[1]
			if qseqid.endswith(basename):
				genes.append(qseqid)
				genomes.append(sseqid.split("|")[-1])   # extract only the genome name

	if not genes:
		logger.warning(f"No matching genes found in BLAST output {blast_output} for genome {basename}.")
		return np.zeros((0,0), dtype=int), [], []

	genes = np.array(genes)
	genomes = np.array(genomes)

	unique_genes, genes_indices = np.unique(genes, return_inverse=True)
	unique_genomes, genomes_indices = np.unique(genomes, return_inverse=True)

	presence_absence = np.zeros((len(unique_genomes), len(unique_genes)), dtype=int)
	np.add.at(presence_absence, (genomes_indices, genes_indices), 1)

	genome_index_map = {genome: i for i, genome in enumerate(unique_genomes)}
	gene_index_map = {gene: j for j, gene in enumerate(unique_genes)}
	return presence_absence, genome_index_map, gene_index_map

def write_pa_table_to_csv(pa_matrix: np.ndarray, genome_map_index: dict[str, int], gene_map_index: dict[str, int], output_path_csv: str) -> None:
	"""
	Write the presence absence matrix to a CSV file with in the genomes as rows and the genes as columns
	Args:
		pa_matrix: NumPy array of shape (num_s_genomes, num_q_genes)
		genome_map_index: A Dictionary mapping the genome name to a row index
		gene_map_index: A Dictionary mapping the gene name to a column index
		output_path_csv: The Path where to save the CSV file
	Returns:
		None
	"""
	if pa_matrix.size == 0:
		logger.warning(f"Empty presence-absence matrix, nothing written to {output_path_csv}.")
		return
	num_genomes, num_genes = pa_matrix.shape
	if len(genome_map_index) != num_genomes or len(gene_map_index) != num_genes:
		logger.error(f"Matrix {pa_matrix} shape does not match mapping.")
		raise ValueError("Matrix {pa_matrix shape doest not match mapping.")
	genomes_sorted = [genome for genome, idx in sorted(genome_map_index.items(), key=lambda x: x[1])]
	genes_sorted = [gene for gene, idx in sorted(gene_map_index.items(), key=lambda x: x[1])]
	try:
		with open(output_path_csv, mode="w", newline='') as file:
			writer = csv.writer(file)
			writer.writerow(["Gene Names"] + genes_sorted)
			for genome in genomes_sorted:
				row_idx = genome_map_index[genome]
				row_data = pa_matrix[row_idx, :]
				writer.writerow([genome] + row_data.tolist())
	except Exception as e:
		logger.error(f"Failed to write presence-absence CSV to {output_path_csv}: {e}.", exc_info=True)
		raise

def calculate_jaccard_similarity_scores(pa_tables: dict[str, dict[np.array, dict, dict]]) -> dict[str, dict[str, float]]:
	"""
	Calculate the Jaccard similarity scores based on a dictioynary containing all the presence-absence tables.
	Args:
		pa_tables: A Dictionary containing all the presence-absence matrices and the corresponding genome and gene mapping dictionaries.
	Returns:
	jaccard_similarity_scores: A Dictionary containing the different Jaccard similarity matrices for each presence-absence table.
	"""
	if not pa_tables:
		logger.error("No presence-absence tables provided for Jaccard calculations.")
		raise ValueError("No presence-absence talbes provided for Jaccard calculations.")
	jaccard_similarity_scores = {}
	for ref_genbank_filename, (pa_matrix, genome_map_index, gene_map_index) in pa_tables.items():
		if pa_matrix.size == 0:
			logger.warning(f"Skipping {ref_genbank_filename}: empty presence-absence matrix.")
			continue
		if ref_genbank_filename not in genome_map_index:
			raise ValueError(f"Reference genome '{ref_genbank_filename}' not found in genome mapping {genome_map_index}.")
		gene_freq = (pa_matrix > 0).sum(axis=0)
		max_freq = gene_freq.max()
		if max_freq == 0:
			logger.warning(f"All gene frequencies are zero for {ref_genbank_filename}, assigning zero weights.")
			weights = np.zeros_like(gene_freq, dtype=float)
		else:
			weights = 1 - (gene_freq / max_freq)
		ref_idx = genome_map_index[ref_genbank_filename]
		ref_vector = pa_matrix[ref_idx, :]
		min_matrix = np.minimum(pa_matrix, ref_vector)
		max_matrix = np.maximum(pa_matrix, ref_vector)
		weighted_intersection = np.sum(weights * min_matrix, axis=1)
		weighted_union = np.sum(weights * max_matrix, axis=1)
		with np.errstate(divide="ignore", invalid="ignore"):
			jaccard_scores = np.where(weighted_union != 0, weighted_intersection / weighted_union, 0.0)
		jac_dict = {genome_name: float(jaccard_scores[idx]) for genome_name, idx in genome_map_index.items()}
		jaccard_similarity_scores[ref_genbank_filename] = jac_dict
	return jaccard_similarity_scores

def construct_matrix(jaccard_similarity_scores: dict[str, dict[str, float]]) -> pd.DataFrame:
	"""
	Construct a complete and symmetrical Jaccard similarity matrix by concatenating and averaging the the Jaccard scores.
	Args:
		jaccard_similarity_scores: A Dictionary containing the different Jaccard similarity matrices for each presence-absence table.
	Returns:
		jaccard_matrix: A Pandas DataFrame containing the complete Jaccard similarity matrix.
	"""
	if not jaccard_similarity_scores:
		logger.error("No Jaccard similarity scores provided.")
		raise ValueError("No Jaccard similarity scores provided.")
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
					logger.warning(f"No similarity scores found for {genome_i} and {genome_j}. Setting to 0.")
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
	try:
		medoid_index = distance_matrix.sum(axis=1).idxmin()
		return medoid_index
	except Exception as e:
		logger.error(f"Failed to compute medoid: {e}")
		raise

def calculate_representation_scores(jaccard_matrix: pd.DataFrame) -> dict[str, float]:
	"""
	Calculate the distances from 1 point to all others which will be the representation scores. The lower the score the better.
	Args:
		jaccard_matrix: A Pandas DataFrame containing the complete Jaccard similarity matrix.
	Returns:
		representation_scores: A Pandas DataFrame containing the representation scores.
	"""
	try:
		distance_matrix = 1 - jaccard_matrix
		rep_scores = distance_matrix.sum(axis=1)
		representation_scores = pd.DataFrame({"Genome": rep_scores.index, "Representation_Score": rep_scores.values})
		representation_scores = representation_scores.sort_values(by="Representation_Score", ascending=True).reset_index(drop=True)
		return representation_scores
	except Exception as e:
		logger.error(f"Failed to calculate representation scores: {e}")
		raise
