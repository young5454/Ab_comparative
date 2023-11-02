import argparse

def fasta_maker_list(faa_file_path, ids_list, save_file_path):
    """
    Makes a new FASTA file by finding ids of ids_list in faa_file
    Saves the new FASTA to save_file_path
    """
    result_sequences = []
    
    # Split by '>' and remove the first empty element
    with open(faa_file_path, 'r') as faa:
        sequences = faa.read().split('>')[1:]  

    count = 0
    faa_ids = []
    
    for seq in sequences:
        ids = seq.split('\n')[0].split(' ')[0]
        faa_ids.append(ids)

        for query in ids_list:
            if query in seq:
                seq_lines = seq.strip().split('\n')
                seq_id = seq_lines[0]

                # Split the fasta entry at the first whitespace
                id_, annotation = seq_id.split(maxsplit=1)

                # Include all
                result_sequences.append(f'>{seq_id}\n')
                result_sequences.append('\n'.join(seq_lines[1:]))
                result_sequences.append('\n')
                count += 1
                
   # Save the result_sequences to a save_file_path
    with open(save_file_path, 'w') as result_file:
        result_file.writelines(result_sequences)
    
    return count


def fasta_maker_text(faa_file_path, ids_text_path, save_file_path):
    """
    Makes a new FASTA file by finding ids of ids_text in faa_file
    Saves the new FASTA to save_file_path
    """
    with open(ids_text_path, 'r') as file:
        ids = file.readlines()
    
    # ids_list: IDs only list
    ids_list = []
    for i in ids:
        woblank = i.strip()
        ids_list.append(woblank)

    result_sequences = []
    
    # Split by '>' and remove the first empty element
    with open(faa_file_path, 'r') as faa:
        sequences = faa.read().split('>')[1:]  

    count = 0
    faa_ids = []
    
    for seq in sequences:
        ids = seq.split('\n')[0].split(' ')[0]
        faa_ids.append(ids)

        for query in ids_list:
            if query in seq:
                seq_lines = seq.strip().split('\n')
                seq_id = seq_lines[0]

                # Split the fasta entry at the first whitespace
                id_, annotation = seq_id.split(maxsplit=1)

                # Include all
                result_sequences.append(f'>{seq_id}\n')
                result_sequences.append('\n'.join(seq_lines[1:]))
                result_sequences.append('\n')
                count += 1
                
   # Save the result_sequences to a save_file_path
    with open(save_file_path, 'w') as result_file:
        result_file.writelines(result_sequences)
    
    return count

def main():    
    # Argument parser
    parser = argparse.ArgumentParser(description='Curate FASTA file based on your need')
    parser.add_argument('--faa_file', required=True, help='Path to Prokka FAA file')
    parser.add_argument('--gene_list', required=True, help='Path to text file of query gene list to be converted')
    parser.add_argument('--save_path', required=False, default='./', help='Path to save log file')
    parser.add_argument('--output', required=False, default='curated.fasta', help='Name of curated FASTA file')
    args = parser.parse_args()

    faa_file = args.faa_file
    gene_list = args.gene_list
    save_path = args.save_path
    output = args.output
    
    path_to_save = save_path + output
    
    with open(path_to_save, 'w') as f:
        fasta_maker_text(faa_file, gene_list, path_to_save)

if __name__ == "__main__":
    main()
