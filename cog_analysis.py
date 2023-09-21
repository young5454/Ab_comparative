import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description='COG Analysis Script')
    parser.add_argument('--tsv_file', required=True, help='Path to annotation TSV file')
    parser.add_argument('--hypo_path', required=True, help='Path to hypothetical proteins fasta file')
    parser.add_argument('--group_name', required=True, help='Name of the strain group to run COG analysis')
    parser.add_argument('--types', required=True, help='Specify core or shell gene types')
    parser.add_argument('--save_path', required=False, help='Path to save pie chart')    
    args = parser.parse_args()

    raw_tsv_file = args.tsv_file
    hypo_path = args.hypo_path
    group_name = args.group_name
    types = args.types
    save_path = args.save_path

    # Count number of hypothetical proteins
    with open(hypo_path, 'r') as hypo:
        sequences = hypo.read().split('>')[1:]  # Split by '>' and remove the first empty element
    
    num_of_hypos = len(sequences)

    # Clean tsv file format for easy parsing
    # Define input and output file names
    cleaned_tsv_file = raw_tsv_file + '.cleaned'

    # Open input and output files
    with open(raw_tsv_file, 'r') as infile, open(cleaned_tsv_file, 'w') as outfile:
        
        # Initialize a flag to skip lines at the beginning
        skip_lines = True
        
        for line in infile:
            # Remove leading and trailing whitespace
            line = line.strip()
       
            # Skip lines that start with "##" at the beginning or are empty
            if line.startswith("##") or not line:
                continue
        
            # Check if the line starts with "#query" and modify it
            if line.startswith("#query"):
                line = line.lstrip("#")
            
            # Write the modified line to the output file
            outfile.write(line + '\n')

    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(cleaned_tsv_file, sep='\t')
    query = df['query']
    cog_cat = df['COG_category']

    # Resolve double, triple categories - weighted
    cog_dictionary_weighted = {}

    # Initialize hypo count
    cog_dictionary_weighted['Hypo'] = num_of_hypos

    # Start counting
    for item in cog_cat:
        if len(item) == 1:
            if item not in cog_dictionary_weighted.keys():
                cog_dictionary_weighted[item] = 1
            elif item in cog_dictionary_weighted.keys():
                curr = cog_dictionary_weighted[item]
                cog_dictionary_weighted[item] = curr + 1
        elif len(item) > 1:
            point = 1 / len(item)
            for bit in item:
                if bit not in cog_dictionary_weighted.keys():
                    cog_dictionary_weighted[bit] = point
                elif bit in cog_dictionary_weighted.keys():
                    curr = cog_dictionary_weighted[bit]
                    cog_dictionary_weighted[bit] = curr + point
    
    # Print result statistics
    print('+--------------------------------------------+')
    print('Number of Hypothetical Proteins are:', num_of_hypos)
    print('+--------------------------------------------+')
    for category in cog_dictionary_weighted.keys():
        if category == 'Hypo':
            continue
        else:
            print('Number of genes in', category + ':', cog_dictionary_weighted[category])
    print('+--------------------------------------------+')
    nums = sum(list(cog_dictionary_weighted.values()))
    print('Total number of queries:', nums)
    
    # Larger grouping: group labels into 4 categories
    grouped = {'Metabolism': 0, 'Information storage and processing': 0, 
            'Mobileome': 0, 'Cellular processing and signaling': 0, 'Hypothetical protein': 0, 'Poorly characterized': 0}

    metabolism = ['C', 'E', 'F', 'G', 'H', 'I', 'P', 'Q']
    info_storage_and_processing = ['A', 'B', 'J', 'K', 'L']
    mobileome = ['X']
    cellular_processing_and_signaling = ['D', 'M', 'N', 'O', 'T', 'U', 'V', 'W', 'Z']
    poorly_characterized = ['R', 'S', '-']

    for ori_key in cog_dictionary_weighted:
        if ori_key in metabolism:
            curr = grouped['Metabolism'] 
            grouped['Metabolism'] = curr + cog_dictionary_weighted[ori_key]
        elif ori_key in info_storage_and_processing:
            curr = grouped['Information storage and processing'] 
            grouped['Information storage and processing'] = curr + cog_dictionary_weighted[ori_key]
        elif ori_key in mobileome:
            curr = grouped['Mobileome']
            grouped['Mobileome'] = curr + cog_dictionary_weighted[ori_key]
        elif ori_key in cellular_processing_and_signaling:
            curr = grouped['Cellular processing and signaling'] 
            grouped['Cellular processing and signaling'] = curr + cog_dictionary_weighted[ori_key]
        elif ori_key == 'Hypo':
            curr = grouped['Hypothetical protein'] 
            grouped['Hypothetical protein'] = curr + cog_dictionary_weighted[ori_key]
        else:
            curr = grouped['Poorly characterized'] 
            grouped['Poorly characterized'] = curr + cog_dictionary_weighted[ori_key]
            
    # Print result statistics
    print('+--------------------------------------------+')
    for category in grouped.keys():
        if category == 'Hypo':
            continue
        else:
            print('Number of genes in', category + ':', grouped[category])
    print('+--------------------------------------------+')
    grouped_nums = sum(list(grouped.values()))
    print('Total number of queries:', grouped_nums) 
    
    print('+--------------------------------------------+')
    print('End of results')
    print('+--------------------------------------------+')

    # Draw pie plot - Grouped
    # Extract labels (categories) and values (proportions) from the dictionary
    labels = list(grouped.keys())
    proportions = list(grouped.values())

    # Define custom colors for the pie chart
    # Metabolism, Information storage and processing, Mobileome, Cellular processing and signaling, Hypo, Poor
    pastel_colors = ['#ff9999', '#ffc000', '#969bff', '#8fd9b6', '#9dd3ff', '#d395d0']             
    explody = [0.03, 0.03, 0.03, 0.03, 0.03, 0.03]

    # Create a pie chart
    plot_title = group_name + ' ' + types + ' ' + 'Genes'
    save_title = group_name + '_' +  types + '_' +  'Genes'
    figure_name = save_title + '.png'
    plt.figure(figsize=(8, 8))
    plt.subplots_adjust(left=0.1)  # You can adjust the value as needed
    plt.pie(proportions, autopct='%1.2f%%', startangle=140, colors=pastel_colors, explode=explody)
    
    plt.title(plot_title)
    plt.axis('equal')  # Equal aspect ratio ensures the pie chart is circular.

    # Add a legend
    plt.legend(bbox_to_anchor=(0.95, 1), loc='upper left', labels=labels)

    # Save pie plot - align to the left with the legends
    figure_name = save_path + figure_name
    plt.savefig(figure_name, bbox_inches='tight')


if __name__ == "__main__":
    main()
