import sys
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
import argparse


def fetch_protein(query):
    # Perform an esearch to get the list of IDs
    handle = Entrez.esearch(db='protein', term=query)
    record = Entrez.read(handle)
    
    # Retrieve the full record for the first ID in the list
    protein_id = record['IdList'][0]
    full_record_handle = Entrez.efetch(db='protein', id=protein_id, rettype='gb', retmode='text')
    full_record = full_record_handle.read()
    
    # Close the handles
    full_record_handle.close()
    
    return full_record
    
    
def name_finder(email, api_key, faa_file, gene_list, n):   
    # NCBI email
    Entrez.email = email

    # NCBI API key
    Entrez.api_key = api_key
    
    # Make query list
    query_list = []
    with open(gene_list, 'r') as query:
        query_list = query.read().split('\n')
    
    # Open faa file and parse data
    with open(faa_file, 'r') as faa:
        # Split by '>' and remove the first empty element
        sequences = faa.read().split('>')[1:]  
    
    for query_id in query_list:
        for sequence in sequences:
            if query_id in sequence:
                seq_lines = sequence.strip().split('\n')
                seq_id = seq_lines[0]  # Extract strain ID & annotation
                sequence = '>' + seq_id + '\n' + '\n'.join(seq_lines[1:])  # Create the modified fasta sequence
                print(sequence)
                print()
    
                # Perform BLASTp on the sequence
                result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
    
                # Parse the BLASTp results
                blast_records = NCBIXML.parse(result_handle)
                cds_symbol_list = []
    
                i = 0  # Initialize the index
                for record in blast_records:
                    if i >= n:
                        break  # Stop when reaching n results
                    for alignment in record.alignments:
                        if i >= n:
                            break  # Stop when reaching n results
                        for hsp in alignment.hsps:
                            i += 1
                            print("+-------------    Alignment #" + str(i) + "    -------------+")
                            alignment_title = alignment.title
                            print("Sequence:", alignment_title)
                            cds_symbol = alignment_title.split('|')[1]
                            cds_symbol_list.append(cds_symbol)
    
                            print("E-value:", hsp.expect)
                            print(hsp.query)
                            print(hsp.match)
                            print(hsp.sbjct)
                            print("+----------    End of Alignment #" + str(i) + "    ---------+")
                            print()
                    print("+--------------------------------------------------+")
                    print("End of BlastP for " + query_id)
                    print("+--------------------------------------------------+")
                    print()

                for cds_symbol in cds_symbol_list:
                    try:
                        pr_full_record = fetch_protein(cds_symbol)
                    except RuntimeError:
                        print('Runtime Error: Skip term')
                        continue
                    except:
                        print('Something else went wrong')
                        continue

                    if '/gene=' in pr_full_record:
                        curr_index = cds_symbol_list.index(cds_symbol) + 1
                        curr_alignment = "#" + str(curr_index)
                        pad = curr_alignment.ljust(5, ' ')
                        print("+--------------------------------------------------+")
                        print(pad + "Protein DB Hit for " + query_id)

                        # Make sure /gene= inside the full record
                        gene_name = pr_full_record.split('/gene=')[1].split()[0].split('"')[1]

                        cds = "CDS symbol: " + cds_symbol
                        gn = "Gene Name: " + str(gene_name)
                        print(cds.rjust(52, ' '))
                        print(gn.rjust(52, ' '))
                        print("+--------------------------------------------------+")
                    else:
                        curr_index = cds_symbol_list.index(cds_symbol) + 1
                        curr_alignment = "#" + str(curr_index)
                        pad = curr_alignment.ljust(5, ' ')
                        print(pad + "No Gene Name found for CDS symbol:", cds_symbol)
                        
                print("+--------------------------------------------------+")
                print("End of results for " + query_id)
                print("+--------------------------------------------------+")
                print()

def main():    
    # Argument parser
    parser = argparse.ArgumentParser(description='Convert Prokka gene ID to NCBI ID considering sequence alignment')
    parser.add_argument('--email', required=True, help='NCBI email is required')
    parser.add_argument('--api_key', required=False, help='NCBI API key if you have one')
    parser.add_argument('--faa_file', required=True, help='Path to Prokka FAA file')
    parser.add_argument('--gene_list', required=True, help='Path to text file of query gene list to be converted')
    parser.add_argument('--n', required=False, default=35, help='Top n number of relevant BlastP alignments. Default 20, Max 50')
    parser.add_argument('--save_path', required=False, default='./', help='Path to save log file')
    parser.add_argument('--output', required=False, default='named.log', help='Name of log file')
    args = parser.parse_args()

    email = args.email
    api_key = args.api_key
    faa_file = args.faa_file
    gene_list = args.gene_list
    n = int(args.n)
    save_path = args.save_path
    output = args.output
    
    path_to_save = save_path + output
    
    with open(path_to_save, 'w') as f:
        # Redirect the standard output to the file
        sys.stdout = f
        name_finder(email, api_key, faa_file, gene_list, n)
        
        # Restore the standard output
        sys.stdout = sys.__stdout__   

    print("+--------------------------------------------------+")
    print("+------    End of results for everything     ------+")
    print("+--------------------------------------------------+")
    
if __name__ == "__main__":
    main()

    
    
