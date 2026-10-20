import argparse
import csv
import re

def parse_dimer_report_to_csv(dimer_report_path, output_csv_path, dg_threshold=None):
	"""
	Parses MFEprimer-3.0 or PrimerPooler dimer report and writes results to CSV.
	CSV headers: primer1, primer2, deltaG, score, Tm, source
	"""
	with open(dimer_report_path, 'r') as f:
		content = f.read()

	results = []
	if "versus" in content:  # PrimerPooler format
		#pattern = r'(\S+)\s+versus\s+(\S+).*?dG\s*=\s*(-?\d+\.\d+)'
		pattern = r'(\S+)\s+versus\s+(\S+).*?dG\s*=\s*(-?\d+\.\d+)\s*$'
		#matches = re.findall(pattern, content)
		matches = re.findall(pattern, content, re.DOTALL | re.MULTILINE)
		for p1, p2, dg in matches:
			dg = float(dg)
			if dg_threshold is None or dg <= dg_threshold:
				results.append({
					"primer1": p1.strip(), 
					"primer2": p2.strip(), 
					"deltaG": dg, 
					"score": "", "Tm": "", 
					"source": "PrimerPooler"
					})
	else:  # MFEprimer format
		pattern = r'Dimer \d+:\s*(.*?)\s*x\s*(.*?)\s+.*?Score:\s*(\d+).*?Tm\s*=\s*(\d+\.\d+).*?Delta G\s*=\s*(-?\d+\.\d+)'
		matches = re.findall(pattern, content)
		for p1, p2, score, tm, dg in matches:
			dg = float(dg)
			if dg_threshold is None or dg <= dg_threshold:
				results.append({
					"primer1": p1.strip(),
					"primer2": p2.strip(),
					"score": int(score),
					"Tm": float(tm),
					"deltaG": dg,
					"source": "MFEprimer"
				})

	# Write to CSV
	headers = ["primer1", "primer2", "deltaG", "score", "Tm", "source"]
	with open(output_csv_path, 'w', newline='') as csvfile:
		writer = csv.DictWriter(csvfile, fieldnames=headers)
		writer.writeheader()
		writer.writerows(results)

	print(f"CSV written to {output_csv_path} with {len(results)} records.")


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Parse dimer report and export to CSV.")
	parser.add_argument("-i", "--input", required=True, help="Path to dimer report file.")
	parser.add_argument("-o", "--output", required=True, help="Path to output CSV file.")
	parser.add_argument("-t", "--threshold", type=float, default=None, help="ΔG threshold (e.g., -5.0).")

	args = parser.parse_args()
	parse_dimer_report_to_csv(args.input, args.output, args.threshold)
 #python parse_dimer_report.py -i dimer_report.txt -o dimer_results.csv -t -5.0
 #python /home2/SardLab/cauret/scripts/gtseq_eval/parse_dimer_report.py -i /home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/mfe_primers_analysis/GTseq_parameters/cisco_100wind_relax_primers_dimer_GTseq_param_mfe.txt -o /home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/dimer_raw_parse/cisco_100wind_relax_primers_dimer_GTseq_param_mfe_parse.csv -t 0
 #python /home2/SardLab/cauret/scripts/gtseq_eval/parse_dimer_report.py -i /home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/GTseq_parameters_pooler/pooler_bonds.txt -o /home2/SardLab/cauret/Ackiss_RADseq_Data/panel_1_files/dimer_raw_parse/pooler_bonds.csv -t 0
 