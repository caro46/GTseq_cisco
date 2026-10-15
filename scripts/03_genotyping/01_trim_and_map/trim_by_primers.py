#!/usr/bin/env python3

"""
Trim FASTA sequences using primer pairs and generate
genomic coordinate annotations.

Input primer table columns
--------------------------
PRIMER_ID
PRIMER_LEFT_SEQUENCE
PRIMER_RIGHT_SEQUENCE

OR

locus
fwd_primer
rev_primer

PRIMER_ID format
----------------
sequence_name_location

Example:
	123_101

where:
	sequence_name = 123
	location = 101

The location is assumed to be the START of the design
window on the catalog sequence and is 1-based.

Stacks FASTA example
--------------------
>1 pos=LG01:1922:+ NS=1
ACTG...

Workflow
--------
1. Match sequence_name to FASTA ID.
2. Convert location from 1-based to 0-based.
3. Search for forward primer near the expected location.
4. Search for reverse primer downstream.
5. Restrict reverse-primer search to the expected
   amplicon size region.
6. Allow mismatches in both primers.
7. Select the best amplicon based on:
	  - mismatch count
	  - expected amplicon size
8. Write trimmed FASTA.
9. Write genomic coordinate table.

Output FASTA headers
--------------------
catalog_id[start:end]

Coordinates in FASTA headers are 1-based inclusive.

Dependencies
------------
pip install biopython pandas

Example
-------
python trim_by_primers.py \
	--primer_table primers.csv \
	--fasta catalog.fa \
	--output trimmed.fa \
	--coordinates trimmed_coordinates.tsv
"""

import argparse
import re
import gzip

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def count_mismatches(seq1, seq2):
	"""
	Count mismatches between two equal-length sequences.
	"""
	return sum(
		a != b
		for a, b in zip(seq1.upper(), seq2.upper())
	)


def find_candidate_matches(
	sequence,
	primer,
	max_mismatches=3,
	start_pos=0,
	end_pos=None
):
	"""
	Find all candidate primer matches.

	Parameters
	----------
	sequence : str
	primer : str
	max_mismatches : int
	start_pos : int
	end_pos : int

	Returns
	-------
	list of tuples

	[
		(position, mismatches),
		...
	]
	"""

	if end_pos is None:
		end_pos = len(sequence)

	primer_len = len(primer)

	matches = []

	for i in range(
		start_pos,
		end_pos - primer_len + 1
	):

		window = sequence[i:i + primer_len]

		mismatches = count_mismatches(
			window,
			primer
		)

		if mismatches <= max_mismatches:
			matches.append(
				(i, mismatches)
			)

	return matches


def parse_stacks_position(description):
	"""
	Parse Stacks FASTA header.

	Example:
		pos=LG01:3439:-

	Returns
	-------
	chromosome, genomic_position, strand
	"""

	match = re.search(
		r"pos=([^:]+):(\d+):([+-])",
		description
	)

	if match is None:
		return None

	chromosome = match.group(1)
	genomic_position = int(match.group(2))
	strand = match.group(3)

	return (
		chromosome,
		genomic_position,
		strand
	)


def main():

	parser = argparse.ArgumentParser(
		description=(
			"Trim FASTA sequences using primer pairs."
		)
	)

	parser.add_argument(
		"--primer_table",
		help=(
			"Primer table (CSV, TSV, or other delimiter-separated text file)"
		),
		required=True
	)

	parser.add_argument(
		"--fasta",
		help="Input FASTA",
		required=True
	)

	parser.add_argument(
		"--output",
		help="Output FASTA",
		required=True
	)

	parser.add_argument(
		"--coordinates",
		default="trimmed_coordinates.tsv",
		help="Coordinate table output"
	)

	parser.add_argument(
		"--max-mismatches",
		type=int,
		default=3,
		help=(
			"Maximum mismatches per primer "
			"(default: 3)"
		)
	)

	parser.add_argument(
		"--location-window",
		type=int,
		default=50,
		help=(
			"Search +/- this many bases around "
			"the PRIMER_ID location "
			"(default: 50)"
		)
	)

	parser.add_argument(
		"--expected-amplicon",
		type=int,
		default=100,
		help=(
			"Expected amplicon size "
			"(default: 100)"
		)
	)

	parser.add_argument(
		"--amplicon-tolerance",
		type=int,
		default=50,
		help=(
			"Allowed deviation from expected "
			"amplicon size "
			"(default: 50)"
		)
	)

	args = parser.parse_args()

	# --------------------------------------------------
	# Read primer table
	# --------------------------------------------------

	try:

		if args.primer_table.lower().endswith(
			".csv"
		):

			primers = pd.read_csv(
				args.primer_table,
				sep=","
			)

		elif args.primer_table.lower().endswith(
			".tsv"
		):

			primers = pd.read_csv(
				args.primer_table,
				sep="\t"
			)

		else:

			primers = pd.read_csv(
				args.primer_table,
				sep=None,
				engine="python"
			)

	except Exception as e:

		raise RuntimeError(
			f"Could not read primer table: {e}"
		)

	primers.columns = (
		primers.columns
		.astype(str)
		.str.strip()
	)

#
# Support multiple primer table formats
#

	# required_columns = {
	# 	"PRIMER_ID",
	# 	"PRIMER_LEFT_SEQUENCE",
	# 	"PRIMER_RIGHT_SEQUENCE"
	# }

	# missing = (
	# 	required_columns
	# 	- set(primers.columns)
	# )

	# if missing:

	# 	raise ValueError(
	# 		"Missing required columns: "
	# 		+ ", ".join(sorted(missing))
	# 	)
	if {
		"PRIMER_ID",
		"PRIMER_LEFT_SEQUENCE",
		"PRIMER_RIGHT_SEQUENCE"
	}.issubset(primers.columns):

		primer_id_col = "PRIMER_ID"
		forward_col = "PRIMER_LEFT_SEQUENCE"
		reverse_col = "PRIMER_RIGHT_SEQUENCE"

	elif {
		"locus",
		"fwd_primer",
		"rev_primer"
	}.issubset(primers.columns):

		primer_id_col = "locus"
		forward_col = "fwd_primer"
		reverse_col = "rev_primer"

	else:

		raise ValueError(
			"Primer table must contain either:\n"
			"  PRIMER_ID, PRIMER_LEFT_SEQUENCE, PRIMER_RIGHT_SEQUENCE\n"
			"or\n"
			"  locus, fwd_primer, rev_primer"
		)

	# --------------------------------------------------
	# Read FASTA (.fa, .fasta, .fa.gz, .fasta.gz)
	# --------------------------------------------------

	# fasta_dict = {
	# 	record.id: record
	# 	for record in SeqIO.parse(
	# 		args.fasta,
	# 		"fasta"
	# 	)
	# }
	if args.fasta.endswith(".gz"):

		fasta_handle = gzip.open(
			args.fasta,
			"rt"
		)

	else:

		fasta_handle = open(
			args.fasta,
			"r"
		)

	fasta_dict = {
		record.id: record
		for record in SeqIO.parse(
			fasta_handle,
			"fasta"
		)
	}

	fasta_handle.close()

	trimmed_records = []

	coordinate_rows = []

	# --------------------------------------------------
	# Process primer pairs
	# --------------------------------------------------

	for _, row in primers.iterrows():

		primer_id = str(
			#row["PRIMER_ID"]
			row[primer_id_col]
		).strip()

		try:

			sequence_name, location = (
				primer_id.rsplit("_", 1)
			)

			#
			# PRIMER_ID location is 1-based
			#
			expected_forward = (
				int(location) - 1
			)

		except ValueError:

			print(
				f"WARNING: Could not parse "
				f"PRIMER_ID '{primer_id}'"
			)

			continue

		if sequence_name not in fasta_dict:

			print(
				f"WARNING: "
				f"{sequence_name} "
				f"not found in FASTA"
			)

			continue

		record = fasta_dict[
			sequence_name
		]

		sequence = str(
			record.seq
		).upper()

		left_primer = str(
			#row["PRIMER_LEFT_SEQUENCE"]
			row[forward_col]
		).upper()

		right_primer = str(
			#row["PRIMER_RIGHT_SEQUENCE"]
			row[reverse_col]
		).upper()

		reverse_primer = str(
			Seq(
				right_primer
			).reverse_complement()
		)

		# ----------------------------------------------
		# Forward primer search window
		# ----------------------------------------------

		forward_window_start = max(
			0,
			expected_forward
			- args.location_window
		)

		forward_window_end = min(
			len(sequence),
			expected_forward
			+ args.location_window
			+ len(left_primer)
		)

		forward_candidates = (
			find_candidate_matches(
				sequence=sequence,
				primer=left_primer,
				max_mismatches=args.max_mismatches,
				start_pos=forward_window_start,
				end_pos=forward_window_end
			)
		)

		if not forward_candidates:

			print(
				f"WARNING: "
				f"No forward primer found for "
				f"{sequence_name}"
			)

			continue

		best_score = float("inf")
		best_interval = None

		# ----------------------------------------------
		# Evaluate candidate amplicons
		# ----------------------------------------------

		for (
			forward_start,
			forward_mismatches
		) in forward_candidates:

			reverse_search_start = (
				forward_start
				+ len(left_primer)
			)

			reverse_search_end = min(
				len(sequence),
				forward_start
				+ args.expected_amplicon
				+ args.amplicon_tolerance
			)

			reverse_candidates = (
				find_candidate_matches(
					sequence=sequence,
					primer=reverse_primer,
					max_mismatches=args.max_mismatches,
					start_pos=reverse_search_start,
					end_pos=reverse_search_end
				)
			)

			if not reverse_candidates:
				continue

			for (
				reverse_start,
				reverse_mismatches
			) in reverse_candidates:

				reverse_end = (
					reverse_start
					+ len(reverse_primer)
				)

				amplicon_length = (
					reverse_end
					- forward_start
				)

				size_penalty = (
					abs(
						amplicon_length
						- args.expected_amplicon
					)
					/
					args.expected_amplicon
				)

				score = (
					forward_mismatches
					+ reverse_mismatches
					+ size_penalty
				)

				if score < best_score:

					best_score = score

					best_interval = (
						forward_start,
						reverse_end
					)

		if best_interval is None:

			print(
				f"WARNING: "
				f"No valid amplicon found for "
				f"{sequence_name}"
			)

			continue

		start_trim, end_trim = best_interval

		trimmed_sequence = sequence[
			start_trim:end_trim
		]

		# ----------------------------------------------
		# Convert Python slice coordinates
		# to 1-based inclusive
		# ----------------------------------------------

		start_1based = start_trim + 1
		end_1based = end_trim

		trimmed_record = record[:]

		trimmed_record.seq = Seq(
			trimmed_sequence
		)

		trimmed_record.id = (
			f"{sequence_name}"
			f"[{start_1based}:{end_1based}]"
		)

		trimmed_record.description = ""

		trimmed_records.append(
			trimmed_record
		)

		# ----------------------------------------------
		# Coordinate output
		# ----------------------------------------------

		parsed = parse_stacks_position(
			record.description
		)

		if parsed is not None:

			chromosome = parsed[0]
			genomic_pos = parsed[1]
			strand = parsed[2]

			if strand == "+":

				genomic_start = (
					genomic_pos
					+ start_trim
				)

				genomic_end = (
					genomic_pos
					+ end_trim
					- 1
				)

			else:

				genomic_start = (
					genomic_pos
					- (end_trim - 1)
				)

				genomic_end = (
					genomic_pos
					- start_trim
				)

			coordinate_rows.append(
				{
					"catalog_id":
						sequence_name,
					"trimmed_id":
						trimmed_record.id,
					"chromosome":
						chromosome,
					"strand":
						strand,
					"genomic_start":
						genomic_start,
					"genomic_end":
						genomic_end,
					"amplicon_length":
						len(trimmed_sequence)
				}
			)

	# --------------------------------------------------
	# Write FASTA
	# --------------------------------------------------

	# SeqIO.write(
	# 	trimmed_records,
	# 	args.output,
	# 	"fasta"
	# )
	if args.output.endswith(".gz"):

		with gzip.open(
			args.output,
			"wt"
		) as out_handle:

			SeqIO.write(
				trimmed_records,
				out_handle,
				"fasta"
			)

	else:

		with open(
			args.output,
			"w"
		) as out_handle:

			SeqIO.write(
				trimmed_records,
				out_handle,
				"fasta"
			)
	# --------------------------------------------------
	# Write coordinate table
	# --------------------------------------------------

	pd.DataFrame(
		coordinate_rows
	).to_csv(
		args.coordinates,
		sep="\t",
		index=False
	)

	print(
		f"Wrote "
		f"{len(trimmed_records)} "
		f"trimmed sequences to "
		f"{args.output}"
	)

	print(
		f"Wrote coordinate table to "
		f"{args.coordinates}"
	)


if __name__ == "__main__":
	main()