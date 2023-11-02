import csv
import argparse

def gene_list_maker(input_csv_file, save_path):
    """
    Makes 6 gene lists for two-strain ROARY results.
    core_all.txt, core_nonhypo.txt, core_hypo.txt
    shells_all.txt, shells_nonhypo.txt, shells_hypo.txt
    """
    
    # Define the input and output file paths
    core_all_file = save_path + 'core_all.txt'
    core_nonhypo_file = save_path + 'core_nonhypo.txt'
    core_hypo_file = save_path + 'core_hypo.txt'

    shells_all_file = save_path + 'shells_all.txt'
    shells_nonhypo_file = save_path + 'shells_nonhypo.txt'
    shells_hypo_file = save_path + 'shells_hypo.txt'

    # Open the CSV file and output files
    with open(input_csv_file, mode='r', newline='') as csv_file, \
            open(core_all_file, mode='w') as core_all, \
            open(core_nonhypo_file, mode='w') as core_nonhypo, \
            open(core_hypo_file, mode='w') as core_hypo, \
            open(shells_all_file, mode='w') as shells_all, \
            open(shells_nonhypo_file, mode='w') as shells_nonhypo, \
            open(shells_hypo_file, mode='w') as shells_hypo:
        csv_reader = csv.reader(csv_file)

        # Skip the header row if present
        next(csv_reader, None)

        for row in csv_reader:
            # Extract the entries starting from column 14
            entries = row[14:]

            # Check if all columns (entries) are filled
            if all(entries):
                core_all.write(row[14] + '\n')
                # Write the leftmost entry (column 14) to core files based on annotation
                if row[2] == 'hypothetical protein':
                    core_hypo.write(entries[0] + '\n')
                else:
                    core_nonhypo.write(entries[0] + '\n')
            else:
                # Select the rightmost entry if shell
                for gene in entries:
                    if gene != '':
                        shell_entry = gene
                shells_all.write(shell_entry + '\n')

                # Check annotation and write to shells files
                if row[2] == 'hypothetical protein':
                    shells_hypo.write(shell_entry + '\n')
                else:
                    shells_nonhypo.write(shell_entry + '\n')
        
        print("Processing complete. Check save_path for results.")

        # for row in csv_reader:
        #     # Check if both columns (15th and 16th, indexed from 0) have data
        #     if row[14] and row[15]:
        #         core_all.write(row[14] + '\n')  # Write to core_all.txt
        #         # Check annotation 
        #         if row[2] == 'hypothetical protein':
        #             core_hypo.write(row[14] + '\n')
        #         else:
        #             core_nonhypo.write(row[14] + '\n')
        #     elif row[14] or row[15]:
        #         shells_all.write(row[14] + row[15] + '\n')  # Write to shells_all.txt
        #         # Check annotation 
        #         if row[2] == 'hypothetical protein':
        #             shells_hypo.write(row[14] + row[15] + '\n')
        #         else:
        #             shells_nonhypo.write(row[14] + row[15] + '\n')
        # print("Processing complete. Check save_path for results.")


def main():    
    # Argument parser
    parser = argparse.ArgumentParser(description='Make CDS core, shells gene list as txt file format')
    parser.add_argument('--csv', required=True, help='Path to gene_presence_absence.csv')
    parser.add_argument('--save_path', required=True, help='Path to save gene lists')
    args = parser.parse_args()

    input_csv_file = args.csv
    save_path = args.save_path
    
    gene_list_maker(input_csv_file, save_path)


if __name__ == "__main__":
    main()

