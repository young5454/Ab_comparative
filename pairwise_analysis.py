import pandas as pd
import argparse

def main():    
    # Argument parser
    parser = argparse.ArgumentParser(description='Type-to-type pairwise comparison analysis')
    parser.add_argument('--csv_file', required=True, help='Path to pairwise blastp output CSV file')
    parser.add_argument('--query_name', required=True, help='Name of the query set')
    parser.add_argument('--subject_name', required=True, help='Name of the subject set')
    parser.add_argument('--evalue', required=False, help='Threshold for evalue')
    parser.add_argument('--pident', required=False, help='Threshold for percent identity')
    parser.add_argument('--bitscore', required=False, help='Threshold for bitscore')
    args = parser.parse_args()
    
    csv_file_path = args.csv_file
    query = args.query_name
    subject = args.subject_name
    evalue_threshold = args.evalue
    pident_threshold = args.pident
    bitscore_threshold = args.bitscore
    
    # Define the new line you want to add
    category = "qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore\n"

    # Read the existing contents of the CSV file
    with open(csv_file_path, 'r') as file:
        first_line = file.readline()
        existing_content = file.read()
        if first_line != category:
            # Concatenate the new line and the existing content
            updated_content = category + existing_content
        else:
            updated_content = existing_content

        # Write the updated content back to the CSV file
        with open(csv_file_path, 'w') as file:
            file.write(updated_content)

#     # Read the TSV file into a pandas DataFrame
#     df = pd.read_csv(csv_file_path, sep='\t')
    
    
if __name__ == "__main__":
    main()
