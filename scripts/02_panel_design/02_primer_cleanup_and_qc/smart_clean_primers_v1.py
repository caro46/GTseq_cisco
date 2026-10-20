import argparse
import csv
import re
from collections import defaultdict

def parse_dimer_report(dimer_report_path, dg_threshold):
	"""
	Parses the MFEprimer-3.0 dimer report and returns:
	- a list of problematic dimers (pairs with ΔG < threshold)
	- a count of how often each primer appears in such dimers
	- a clean set of primers
	"""
	problematic_dimers = []
	problematic_dimers_dg = []
	primer_dimer_count = defaultdict(int)
	#defaultdict automatically creates a default value for any key that doesn’t exist yet.

	# Read the entire dimer report file
	with open(dimer_report_path, 'r') as f:
		content = f.read()

	# Use regex to extract primer pairs and their ΔG values
	# Format expected: 
	# - MFEprimer: "Dimer #: primer1 x primer2 ... Delta G = -X.XX"
	# - primerpooler "primer1 versus primer2 ... dG = -X.XX"

	# Determine which pattern to use based on file content or metadata

	if "versus" in content:  #primerpooler
		#pattern = r'(.*?) versus (.*?)\s+.*?dG = (-?\d+\.\d+).*?(?=\n\n|\Z)'
		pattern = r'(\S+)\s+versus\s+(\S+).*?dG\s*=\s*(-?\d+\.\d+)\s*$'
	else:
		pattern = r'Dimer \d+: (.*?) x (.*?)\s+.*?Delta G = (-?\d+\.\d+)'

	#searches a string for all non-overlapping matches of a pattern and returns them as a list.
	#dimer_matches = re.findall(pattern, content, re.DOTALL)
	dimer_matches = re.findall(pattern, content, re.DOTALL | re.MULTILINE)
	#return matches

	#dimer_matches = re.findall( 
	#	#searches a string for all non-overlapping matches of a pattern and returns them as a list.
	#	r'Dimer \d+: (.*?) x (.*?)\s+.*?Delta G = (-?\d+\.\d+)', content, re.DOTALL
	#)

	for primer1, primer2, dg in dimer_matches:
		dg = float(dg)
		if dg < dg_threshold:
			# If ΔG is below threshold, consider this a problematic dimer
			problematic_dimers.append((primer1.strip(), primer2.strip()))
			problematic_dimers_dg.append((primer1.strip(), primer2.strip(),dg))
			# Count how often each primer appears in problematic dimers
			primer_dimer_count[primer1.strip()] += 1
			primer_dimer_count[primer2.strip()] += 1

	return problematic_dimers, primer_dimer_count, problematic_dimers_dg

def iterative_primer_removal(problematic_dimers):
	"""
	Iteratively removes the most problematic primer (most frequent in dimers),
	and removes all dimers involving it, until no problematic dimers remain.
	Returns a set of primers to remove.
	"""
	removed_primers = set()
	dimers = problematic_dimers.copy()

	while dimers:
		# Count how many times each primer appears in current dimer list
		primer_count = defaultdict(int)
		for p1, p2 in dimers:
			primer_count[p1] += 1
			primer_count[p2] += 1

		# Find the primer involved in the most dimers
		most_problematic = max(primer_count.items(), key=lambda x: x[1])[0]
		removed_primers.add(most_problematic)

		# Remove all dimers that include the most problematic primer
		dimers = [(p1, p2) for p1, p2 in dimers if most_problematic not in (p1, p2)]

	return removed_primers

def filter_primer_list(primer_list_path, removed_primers, output_path):
	"""
	Filters the input primer list file (tab-separated marker name, forward/reverse columns),
	and writes only those pairs where both primers are retained.
	"""
	with open(primer_list_path, 'r') as infile, open(output_path, 'w', newline='') as outfile:
		reader = csv.reader(infile, delimiter='\t')
		writer = csv.writer(outfile, delimiter='\t')

		for row in reader:
			if len(row) < 3:
				continue  # Skip malformed lines
			marker, fwd, rev = row[0].strip(), row[1].strip(), row[2].strip() #0 is marker name
			# Only write pairs where neither primer was removed
			if fwd not in removed_primers and rev not in removed_primers:
				writer.writerow([marker,fwd, rev])

#def write_problematic_dimers(problematic_dimers, output_path):
def write_problematic_dimers(problematic_dimers_dg, output_path):
	"""
	Writes the list of problematic dimers to a CSV file.
	"""
	with open(output_path, 'w', newline='') as f:
		writer = csv.writer(f)
		writer.writerow(['Primer1', 'Primer2', 'DeltaG'])
		for p1, p2, dg in problematic_dimers_dg:
			writer.writerow([p1, p2, dg])

def write_primer_dimer_count(primer_dimer_count, output_path):
	"""
	Writes the primer dimer count dictionary to a CSV file.
	"""
	with open(output_path, 'w', newline='') as f:
		writer = csv.writer(f)
		writer.writerow(['Primer', 'DimerCount'])
		for primer, count in primer_dimer_count.items():
			writer.writerow([primer, count])


def main():
	# Set up command-line argument parsing
	parser = argparse.ArgumentParser(description='Filter primers based on MFEprimer dimer report.')
	parser.add_argument('--dimer_report', required=True, help='Path to MFEprimer dimer report file')
	parser.add_argument('--primer_list', required=True, help='Path to primer list file (tab-separated)')
	parser.add_argument('--output', required=True, help='Path to save cleaned primer list')
	parser.add_argument('--dg_threshold', type=float, default=-9.0, help='Delta G threshold (default: -9.0 kcal/mol)')
	parser.add_argument('--problematic_output', default='problematic_dimers.csv', help='Output file for problematic dimers')
	parser.add_argument('--count_output', default='primer_dimer_counts.csv', help='Output file for primer dimer counts')

	args = parser.parse_args()

	# Step 1: Parse dimer report to find problematic dimers and primer frequencies
	problematic_dimers, primer_dimer_count, problematic_dimers_dg = parse_dimer_report(args.dimer_report, args.dg_threshold)

	# Step 2: Iteratively remove primers most involved in problematic dimers
	removed_primers = iterative_primer_removal(problematic_dimers)

	# Step 3: Filter the original primer list to exclude removed primers
	filter_primer_list(args.primer_list, removed_primers, args.output)

	# Step 4: Write outputs for problematic dimers and count of dimers for each primer
	write_problematic_dimers(problematic_dimers_dg, args.problematic_output)
	write_primer_dimer_count(primer_dimer_count, args.count_output)

	# Step 4: Print summary of results
	print(f"✅ Cleaned primer list saved to {args.output}")
	print(f"🔍 Removed {len(removed_primers)} problematic primers.")
	if removed_primers:
		# Identify the single most problematic primer from the original count
		most_problematic = max(primer_dimer_count.items(), key=lambda x: x[1])
		print(f"🔥 Most problematic primer: {most_problematic[0]} (involved in {most_problematic[1]} dimers)")

# Example usage:
# python smart_clean_primers.py \
#   --dimer_report dimer_report.txt \
#   --primer_list all_primers.tsv \
#   --output cleaned_primers.tsv \
#   --dg_threshold -9.0

# python /home2/SardLab/cauret/scripts/gtseq_eval/smart_clean_primers.py \
# --dimer_report /home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/mfe_primers_analysis/cisco_100wind_relax_primers_dimer_mfe.txt \
# --primer_list /home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/marker_real_primer_set.txt \
# --output /home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/clean_primer_set/clean_primers.txt \
# --problematic_output /home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/clean_primer_set/problematic_dimers.csv \
# --count_output /home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/clean_primer_set/primer_dimer_counts.csv

if __name__ == '__main__':
	main()
