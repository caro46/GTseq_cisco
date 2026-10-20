import re
import argparse

def extract_problematic_primers(dimer_report_path, dg_threshold):
    """
    Parses the MFEprimer-3.0 dimer report and returns a set of primers involved in dimers
    with Delta G below the specified threshold.
    """
    problematic_primers = set()
    with open(dimer_report_path, 'r') as f:
        content = f.read()
        dimer_matches = re.findall(
#            r'Dimer \\d+: (.*?) x (.*?)\\s+.*?Delta G = (-?\\d+\\.\\d+)', content, re.DOTALL
            r'Dimer \d+: (.*?) x (.*?)\s+.*?Delta G = (-?\d+\.\d+)', content, re.DOTALL

        )
        for primer1, primer2, dg in dimer_matches:
            if float(dg) < dg_threshold:
                problematic_primers.update([primer1.strip(), primer2.strip()])
    return problematic_primers

def load_all_primers(primer_list_path):
    """
    Loads all primers from a file. Assumes one primer ID per line.
    """
    with open(primer_list_path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def save_cleaned_primers(cleaned_primers, output_path):
    """
    Saves the cleaned list of primers to a file.
    """
    with open(output_path, 'w') as f:
        for primer in cleaned_primers:
            f.write(primer + '\\n')

def main():
    parser = argparse.ArgumentParser(description='Clean primer list based on MFEprimer dimer report.')
    parser.add_argument('--dimer_report', required=True, help='Path to MFEprimer dimer report file')
    parser.add_argument('--primer_list', required=True, help='Path to file containing all primer IDs')
    parser.add_argument('--output', required=True, help='Path to save cleaned primer list')
    parser.add_argument('--dg_threshold', type=float, default=-9.0, help='Delta G threshold (default: -9.0 kcal/mol)')

    args = parser.parse_args()

    problematic_primers = extract_problematic_primers(args.dimer_report, args.dg_threshold)
    all_primers = load_all_primers(args.primer_list)
    cleaned_primers = [p for p in all_primers if p not in problematic_primers]
    save_cleaned_primers(cleaned_primers, args.output)

    print(f"✅ Cleaned primer list saved to {args.output}.")
    print(f"🔍 {len(cleaned_primers)} primers retained out of {len(all_primers)}.")

#ex to run
#python clean_primers.py \
#  --dimer_report dimer_report.txt \
#  --primer_list all_primers.txt \
#  --output cleaned_primers.txt \
#  --dg_threshold -6.0

if __name__ == '__main__':
    main()
