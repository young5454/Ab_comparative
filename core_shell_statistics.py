import argparse
from fasta_maker import fasta_maker_text

def core_all_fasta(faa_file, text_path, save_path):
    """
    Curates core_all.fasta from core_all.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """

    query_id = text_path + 'core_all.txt'
    save_path = save_path + 'core_all.fasta'

    with open(query_id, 'r') as q:
        query_ids = q.readlines()
        
    # Extract IDs only
    final_query = []
    for query_id in query_ids:
        woblank = query_id.strip()
        final_query.append(woblank)

    result_sequences = []
    no_entry_cores = []

    with open(faa_file, 'r') as faa:
        sequences = faa.read().split('>')[1:]  # Split by '>' and remove the first empty element

    faa_ids = []
    cds_all_core_count = 0

    for seq in sequences:
        ids = seq.split('\n')[0].split(' ')[0]
        faa_ids.append(ids)
        for query in final_query:
            if query in seq:
                seq_lines = seq.strip().split('\n')
                seq_id = seq_lines[0]
                
                # Split the fasta entry at the first whitespace
                id_, annotation = seq_id.split(maxsplit=1)
                
                # Include all
                result_sequences.append(f'>{seq_id}\n')
                result_sequences.append('\n'.join(seq_lines[1:]))
                result_sequences.append('\n') 
                cds_all_core_count += 1
    
    # These are core genes that are not protein-coding
    non_cds_core_count = 0
    for query in final_query:
        if query not in faa_ids:
            no_entry_cores.append(query)
            non_cds_core_count += 1

    # Save the result_sequences to a file
    with open(save_path, 'w') as result_file:
        result_file.writelines(result_sequences)
    
    return cds_all_core_count, non_cds_core_count, no_entry_cores


def core_nonhypo_fasta(faa_file, text_path, save_path):
    """
    Curates core_nonhypo.fasta from core_nonhypo.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """

    query_id = text_path + 'core_nonhypo.txt'
    save_path = save_path + 'core_nonhypo.fasta'

    cds_core_non_hypo_count = fasta_maker_text(faa_file, query_id, save_path)

    return cds_core_non_hypo_count


def core_hypo_fasta(faa_file, text_path, save_path):
    """
    Curates core_hypo.fasta from core_hypo.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """

    query_id = text_path + 'core_hypo.txt'
    save_path = save_path + 'core_hypo.fasta'

    cds_core_hypo_count = fasta_maker_text(faa_file, query_id, save_path)
    
    return cds_core_hypo_count


def shells_all_fasta(all_faa_files, text_path, save_path):
    """
    Curates shells_all.fasta from shells_all.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """

    query_id = text_path + 'shells_all.txt'
    save_path = save_path + 'shells_all.fasta'

    with open(query_id, 'r') as q:
        query_ids = q.readlines()
        
    # Extract IDs only
    final_query = []
    for query_id in query_ids:
        woblank = query_id.strip()
        final_query.append(woblank)

    result_sequences = []
    no_entry_shells = []

    cds_all_shells_count = 0
    non_cds_shells_count = 0

    for faa_file in all_faa_files:
        with open(faa_file, 'r') as faa:
            sequences = faa.read().split('>')[1:]  # Split by '>' and remove the first empty element
        
        faa_ids = [] 
        
        for seq in sequences:
            ids = seq.split('\n')[0].split(' ')[0]
            faa_ids.append(ids)
            
            for query in final_query:
                if query in seq:
                    seq_lines = seq.strip().split('\n')
                    seq_id = seq_lines[0]
                
                    # Split the fasta entry at the first whitespace
                    id_, annotation = seq_id.split(maxsplit=1)
                    result_sequences.append(f'>{seq_id}\n')
                    result_sequences.append('\n'.join(seq_lines[1:]))
                    result_sequences.append('\n') 
                    cds_all_shells_count += 1

        # These are shell genes that are not protein-coding
        curr_strain = faa_ids[0].split('_')[1]
        for query in final_query :
            strain = query.split('_')[1]
            if query not in faa_ids:
                if strain == curr_strain:
                    no_entry_shells.append(query)
                    non_cds_shells_count += 1


    # Save the result_sequences to a file
    with open(save_path, 'w') as result_file:
        result_file.writelines(result_sequences)
    
    return cds_all_shells_count, non_cds_shells_count, no_entry_shells
    

def shells_nonhypo_fasta(all_faa_files, text_path, save_path):
    """
    Curates shells_nonhypo.fasta from shells_nonhypo.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """
    query_id = text_path + 'shells_nonhypo.txt'
    save_path = save_path + 'shells_nonhypo.fasta'

    with open(query_id, 'r') as q:
        query_ids = q.readlines()
        
    # Extract IDs only
    final_query = []
    for query_id in query_ids:
        woblank = query_id.strip()
        final_query.append(woblank)

    result_sequences = []

    cds_shells_non_hypo_count = 0

    for faa_file in all_faa_files:
        with open(faa_file, 'r') as faa:
            sequences = faa.read().split('>')[1:]  # Split by '>' and remove the first empty element
        
        faa_ids = [] 
        
        for seq in sequences:
            ids = seq.split('\n')[0].split(' ')[0]
            faa_ids.append(ids)
            
            for query in final_query:
                if query in seq:
                    seq_lines = seq.strip().split('\n')
                    seq_id = seq_lines[0]
                
                    # Split the fasta entry at the first whitespace
                    id_, annotation = seq_id.split(maxsplit=1)
                    result_sequences.append(f'>{seq_id}\n')
                    result_sequences.append('\n'.join(seq_lines[1:]))
                    result_sequences.append('\n')
                    cds_shells_non_hypo_count += 1

    # Save the result_sequences to a file
    with open(save_path, 'w') as result_file:
        result_file.writelines(result_sequences)
    
    return cds_shells_non_hypo_count


def shells_hypo_fasta(all_faa_files, text_path, save_path):
    """
    Curates shells_hypo.fasta from shells_hypo.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """
    query_id = text_path + 'shells_hypo.txt'
    save_path = save_path + 'shells_hypo.fasta'
    with open(query_id, 'r') as q:
        query_ids = q.readlines()
        
    # Extract IDs only
    final_query = []
    for query_id in query_ids:
        woblank = query_id.strip()
        final_query.append(woblank)

    result_sequences = []

    cds_shells_hypo_count = 0

    for faa_file in all_faa_files:
        with open(faa_file, 'r') as faa:
            sequences = faa.read().split('>')[1:]  # Split by '>' and remove the first empty element
        
        faa_ids = [] 
        
        for seq in sequences:
            ids = seq.split('\n')[0].split(' ')[0]
            faa_ids.append(ids)
            
            for query in final_query:
                if query in seq:
                    seq_lines = seq.strip().split('\n')
                    seq_id = seq_lines[0]
                
                    # Split the fasta entry at the first whitespace
                    id_, annotation = seq_id.split(maxsplit=1)
                    result_sequences.append(f'>{seq_id}\n')
                    result_sequences.append('\n'.join(seq_lines[1:]))
                    result_sequences.append('\n')
                    cds_shells_hypo_count += 1

    # Save the result_sequences to a file
    with open(save_path, 'w') as result_file:
        result_file.writelines(result_sequences)
    
    return cds_shells_hypo_count


def main():    
    # Argument parser
    parser = argparse.ArgumentParser(description='Curate FASTAs of CDS core, shells genes. Also print result statistics')
    parser.add_argument('--text_path', required=True, help='Path to gene list. Never omit a / at the end')
    parser.add_argument('--save_path', required=True, help='Path to save gene FASTAs. Never omit a / at the end')
    parser.add_argument('--all_faas', required=True, help='List of all faa file paths in current group. Format: path1$path2$path3.. in one single string')
    parser.add_argument('--all_roary', required=False, help='Number of all ROARY genes, annotated from PROKKA')
    parser.add_argument('--all_roary_core', required=True, help='Number of all ROARY core genes, including non-protein-coding genes')
    parser.add_argument('--all_roary_shell', required=True, help='Number of all ROARY shell genes, including non-protein-coding genes')
    args = parser.parse_args()

    text_path = args.text_path
    save_path = args.save_path
    all_faas = args.all_faas
    all_roary_genes = args.all_roary
    roary_core_genes = int(args.all_roary_core)
    roary_shell_genes = int(args.all_roary_shell)

    # Split input faa_files into a list type
    all_faa_files = all_faas.split('$')
    faa_file = all_faa_files[0]

    # Function calls
    cds_all_core_count, non_cds_core_count, no_entry_cores = core_all_fasta(faa_file, text_path, save_path)
    cds_core_non_hypo_count = core_nonhypo_fasta(faa_file, text_path, save_path)
    cds_core_hypo_count = core_hypo_fasta(faa_file, text_path, save_path)
    
    cds_all_shells_count, non_cds_shells_count, no_entry_shells = shells_all_fasta(all_faa_files, text_path, save_path)
    cds_shells_non_hypo_count = shells_nonhypo_fasta(all_faa_files, text_path, save_path)
    cds_shells_hypo_count = shells_hypo_fasta(all_faa_files, text_path, save_path)
    
    # Print Result Statistics
    print('+--------------------------------------------------------------------+')
    print('+                          RESULT STATISTICS                         +')
    print('+--------------------------------------------------------------------+')
    print('+--------------------------------------------------------------------+')
    print('+        Successfully wrote all core queries to result FASTAs        +')
    print('+--------------------------------------------------------------------+')
    print('The total number of ROARY core genes are:', roary_core_genes)
    print('>> Please refer to summary_statistics.txt for this info \n')
    print('The total number of CDS all core genes are:', cds_all_core_count)
    print('The total number of CDS non-hypo core genes are:', cds_core_non_hypo_count)
    print('The total number of CDS hypo core genes are:', cds_core_hypo_count)
    print('The total number of non-CDS core genes are:', non_cds_core_count)
    print('+--------------------------------------------------------------------+')

    print('# of CDS non-hypo core + # of CDS hypo core = # of CDS all core?')
    sum = cds_core_non_hypo_count + cds_core_hypo_count
    print('>>', sum == cds_all_core_count)
    print('+--------------------------------------------------------------------+')

    print('# of CDS all core + # of non-CDS core = # of ROARY core?')
    sumy = cds_all_core_count + non_cds_core_count
    print('>>', sumy == roary_core_genes)
    print('+--------------------------------------------------------------------+')

    print('The ids of non-CDS core genes are as follows:')
    for ids in no_entry_cores: print(ids)
    print('These are one of the followings: tRNA, rRNA, tmRNA, ncRNA, ...')
    print('>> Please refer to .ffn file for each gene info')
        
    print('+--------------------------------------------------------------------+')
    print('+       Successfully wrote all shell queries to result FASTAs        +')
    print('+--------------------------------------------------------------------+')
    print('The total number of ROARY shell genes are:', roary_shell_genes)
    print('>> Please refer to summary_statistics.txt for this info \n')
    print('The total number of CDS all shell genes are:', cds_all_shells_count)
    print('The total number of CDS non-hypo shell genes are:', cds_shells_non_hypo_count)
    print('The total number of CDS hypo shell genes are:', cds_shells_hypo_count)
    print('The total number of non-CDS shell genes are:', non_cds_shells_count)
    print('+--------------------------------------------------------------------+')

    print('# of CDS non-hypo shells + # of CDS hypo shells = # of CDS all shells?')
    sumyy = cds_shells_non_hypo_count + cds_shells_hypo_count
    print('>>', sumyy == cds_all_shells_count)
    print('+--------------------------------------------------------------------+')

    print('# of CDS all shells + # of non-CDS shells = # of ROARY shells?')
    sumyyy = cds_all_shells_count + non_cds_shells_count
    print('>>', sumyyy == roary_shell_genes)
    print('+--------------------------------------------------------------------+')

    print('The ids of non-CDS shell genes are as follows:')
    for ids in no_entry_shells: print(ids)
    print('These are one of the followings: tRNA, rRNA, tmRNA, ncRNA, ...')
    print('>> Please refer to .ffn file for each gene info')
    print('+--------------------------------------------------------------------+')


if __name__ == "__main__":
    main()