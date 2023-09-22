import pandas as pd
import argparse
from fasta_maker import fasta_maker_list
    
def main():    
    # Argument parser
    parser = argparse.ArgumentParser(description='Type-to-type protein pairwise comparison analysis')
    parser.add_argument('--csv_file', required=True, help='Path to pairwise blastp output CSV file')
    parser.add_argument('--query_name', required=True, help='Name of the query protein set')
    parser.add_argument('--subject_name', required=True, help='Name of the subject protein set')
    parser.add_argument('--query_faa_file', required=True, help='Path to representative faa file of query proteins')
    parser.add_argument('--subject_faa_file', required=True, help='Path to representative faa file of subject proteins')
    # parser.add_argument('--analysis_type', required=True, choices=['core-to-core', 'core-to-strain'], \
                        # help='Choose between core-to-core or core-to-strain')
    parser.add_argument('--evalue', required=False, default=80, help='Threshold for evalue. Default value is 80')
    parser.add_argument('--pident', required=False, default=0.001, help='Threshold for percent identity. Default value is 0.001')
    parser.add_argument('--save_query', required=False, default='./', help='Path to save query pass FASTA. Defualt is current directory')
    parser.add_argument('--save_subject', required=False, default='./', help='Path to save subject pass FASTA. Defualt is current directory')
    parser.add_argument('--output_query', required=False, default='pass_query.fasta', help='File name of query pass FASTA. Defualt is pass_query.fasta')
    parser.add_argument('--output_subject', required=False, default='pass_subject.fasta', help='File name of subject pass FASTA. Defualt is pass_subject.fasta')
    args = parser.parse_args()
    
    csv_file_path = args.csv_file
    
    query = args.query_name
    subject = args.subject_name
    
    query_faa_file = args.query_faa_file
    subject_faa_file = args.subject_faa_file    
    
    # analysis_type = args.analysis_type
    evalue_threshold = float(args.evalue)
    pident_threshold = float(args.pident)
    # bitscore_threshold = args.bitscore
    
    save_query = args.save_query
    save_subject = args.save_subject
    
    output_query = args.output_query
    output_subject = args.output_subject
    
    # Concatenate categories to the first line of CSV file
    category = "qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore\n"

    # Read the existing contents of the CSV file
    with open(csv_file_path, 'r') as file:
        first_line = file.readline()
        existing_content = first_line + file.read()
        
        if first_line == category:
            updated_content = existing_content
        else:
            updated_content = category + existing_content

        # Write the updated content back to the CSV file
        with open(csv_file_path, 'w') as file:
            file.write(updated_content)

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_file_path, sep='\t')
    query = df['qseqid']
    subject = df['sseqid']
    evalue = df['evalue']
    identity = df['pident']
    
    # Count query/subject entries that pass the identity threshold & eval threshold
    # Currently support only percent identity and e-value
    above_threshold = []
    pass_queries = []
    pass_subjects = []
    
    for i in range(len(identity)):
        curr_identity = identity[i]
        curr_evalue = evalue[i]
        
        if (curr_identity >= pident_threshold) and (curr_evalue < evalue_threshold):
            pass_query = query[i]
            pass_subject = subject[i]
            pass_queries.append(pass_query)
            pass_subjects.append(pass_subject)

    # Remove redundant entries
    pass_query_set = set(pass_queries)
    pass_subject_set = set(pass_subjects)

    print('The number of unique query entries that pass the set threshold:', len(pass_query_set))
    print('The number of unique subject entries that pass the set threshold:', len(pass_subject_set))
    
    # Make FASTAs of passed entries (pass FASTAs)
    # Pass FASTA - query
    query_save_file_name = save_query + output_query
    query_count = fasta_maker_list(query_faa_file, pass_query_set, query_save_file_name)
    print('Successfully created and saved query pass FASTA file')
    
    # Pass FASTA - subject
    subject_save_file_name = save_subject + output_subject
    subject_count = fasta_maker_list(subject_faa_file, pass_subject_set, subject_save_file_name) 
    print('Successfully created and saved subject pass FASTA file')

    
if __name__ == "__main__":
    main()
