import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import argparse
import csv
import numpy as np

def main():    
    # Argument parser
    parser = argparse.ArgumentParser(description='COG Analysis Script')
    parser.add_argument('--tsv_file', required=True, help='Path to annotation TSV file')
    parser.add_argument('--hypo_path', required=True, help='Path to hypothetical proteins FASTA file')
    parser.add_argument('--group_name', required=True, help='Name of the strain group to run COG analysis')
    parser.add_argument('--types', required=True, help='Specify core or shell gene types')
    parser.add_argument('--save_path', required=False, help='Path to save pie chart')
    parser.add_argument('--num_labels', required=False, default=False, help='Show weighted counts for COG categories')    
    args = parser.parse_args()

    # Define argument variables
    raw_tsv_file = args.tsv_file
    hypo_path = args.hypo_path
    group_name = args.group_name
    types = args.types
    save_path = args.save_path
    num_labels = args.num_labels

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

    # Read the cleaned TSV file into a pandas DataFrame
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
    print('Number of Hypothetical proteins are:', num_of_hypos)
    print('+--------------------------------------------+')

    for category in cog_dictionary_weighted.keys():
        if category == 'Hypo':
            continue
        else:
            # Print weighted count in three decimal places
            print('Number of genes in', category + ':', "{:.3f}".format(cog_dictionary_weighted[category]))
    print('+--------------------------------------------+')
    nums = sum(list(cog_dictionary_weighted.values()))
    print('Total number of queries:', "{:.1f}".format(nums))

    # Larger grouping: group labels into 4 categories
    grouped = {'Cellular processing and signaling': 0, 'Hypothetical protein': 0,  'Information storage and processing': 0, 'Metabolism': 0,
            'Mobileome': 0, 'Poorly characterized': 0}

    metabolism = ['C', 'E', 'F', 'G', 'H', 'I', 'P', 'Q']
    info_storage_and_processing = ['A', 'B', 'J', 'K', 'L']
    mobileome = ['X']
    cellular_processing_and_signaling = ['D', 'M', 'N', 'O', 'T', 'U', 'V', 'W', 'Z']
    poorly_characterized = ['R', 'S', '-']
    hypo = ['Hypo']

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
            print('Number of genes in', category + ':', "{:.3f}".format(grouped[category]))
    print('+--------------------------------------------+')
    grouped_nums = sum(list(grouped.values()))
    print('Total number of queries:', "{:.1f}".format(grouped_nums))
    print('+--------------------------------------------+')
    print('End of results')
    print('+--------------------------------------------+')

    # Create COG CSV file
    cogfile_name = group_name + '_' + types + '_' + 'cog.csv'
    with open(cogfile_name, 'w') as cogfile:
        cogfile.write('group,cog,counts\n')
        sequence = [cellular_processing_and_signaling, hypo, info_storage_and_processing, metabolism, mobileome, poorly_characterized]
        for item in cellular_processing_and_signaling:
            if item in cog_dictionary_weighted.keys():
                weighted_count = cog_dictionary_weighted[item]
                cogfile.write('Cellular processing and signaling,' + str(item) + ',' + str(weighted_count) + '\n')
            else:
                cogfile.write('Cellular processing and signaling,' + str(item) + ','  + '0' + '\n')

        for item in hypo:
            if item in cog_dictionary_weighted.keys():
                weighted_count = cog_dictionary_weighted[item]
                cogfile.write('Hypothetical protein,' + str(item) + ',' + str(weighted_count) + '\n')
            else:
                cogfile.write('Hypothetical protein,' + str(item) + ','  + '0' + '\n')
        
        for item in info_storage_and_processing:
            if item in cog_dictionary_weighted.keys():
                weighted_count = cog_dictionary_weighted[item]
                cogfile.write('Information storage and processing,' + str(item) + ',' + str(weighted_count) + '\n')
            else:
                cogfile.write('Information storage and processing,' + str(item) + ','  + '0' + '\n')
        
        for item in metabolism:
            if item in cog_dictionary_weighted.keys():
                weighted_count = cog_dictionary_weighted[item]
                cogfile.write('Metabolism,' + str(item) + ',' + str(weighted_count) + '\n')
            else:
                cogfile.write('Metabolism,' + str(item) + ','  + '0' + '\n')
        
        for item in mobileome:
            if item in cog_dictionary_weighted.keys():
                weighted_count = cog_dictionary_weighted[item]
                cogfile.write('Mobileome,' + str(item) + ',' + str(weighted_count) + '\n')
            else:
                cogfile.write('Mobileome,' + str(item) + ','  + '0' + '\n')

        for item in poorly_characterized:
            if item in cog_dictionary_weighted.keys():
                weighted_count = cog_dictionary_weighted[item]
                cogfile.write('Poorly characterized,' + str(item) + ',' + str(weighted_count) + '\n')
            else:
                cogfile.write('Poorly characterized,' + str(item) + ','  + '0' + '\n')
        

    # Plot nested pie chart of COG results
    df = pd.read_csv(cogfile_name)
    vals = df['counts']
    groupsum = df.groupby('group')['counts'].sum()

    facecolor = '#ffffff'
    font_color = '#000000'
    groups_label = ['Cellular processing \n and signaling', 'Hypothetical protein',
                    'Information storage \n and processing', 'Metabolism',  \
                    'Mobileome',
                    'Poorly characterized']

    cogs_label = ['D', 'M', 'N', 'O', 'T', 'U', 'V', 'W', 'Z',  # Cellular processing and signaling
                  'Hypo',                                       # Hypothetical protein
                  'A', 'B', 'J', 'K', 'L',                      # Information storage and processing
                  'C', 'E', 'F', 'G', 'H', 'I', 'P', 'Q',       # Metabolism
                  'X',                                          # Mobileome
                  'R', 'S', '-']                                # Poorly characterized

    # Only label those that have nonzero values
    cogs_label_nonzero = []
    cogs_label_nonzero_num = []
    for label in cogs_label:
        if label in cog_dictionary_weighted.keys():
            cogs_label_nonzero.append(label)
            count = cog_dictionary_weighted[label]
            label = label + ' ' + '(' + "{:.3f}".format(count) + ')'
            cogs_label_nonzero_num.append(label)
        else:
            cogs_label_nonzero.append('')
            cogs_label_nonzero_num.append('')

    groups_legend = ['Cellular processing and signaling', 'Hypothetical protein',
                    'Information storage and processing', 'Metabolism',  \
                    'Mobileome',
                    'Poorly characterized']

    cogs_legend = ['[D] Cell cycle control, cell division, chromosome partitioning',
                   '[M] Cell wall/membrane/envelope biogenesis',
                   '[N] Cell motility', 
                   '[O] Posttranslational modification, protein turnover, chaperones', 
                   '[T] Signal transduction mechanisms', 
                   '[U] Intracellular trafficking, secretion, and vesicular transport', 
                   '[V] Defense mechanisms', 
                   '[W] Extracellular structures', 
                   '[Z] Cytoskeleton',
                   '[Hypo] Hypothetical proteins',
                   '[A] RNA processing and modification', 
                   '[B] Chromatin structure and dynamics', 
                   '[J] Translation, ribosomal structure and biogenesis', 
                   '[K] Transcription', 
                   '[L] Replication, recombination and repair',
                   '[C] Energy production and conversion', 
                   '[E] Amino acid transport and metabolism', 
                   '[F] Nucleotide transport and metabolism', 
                   '[G] Carbohydrate transport and metabolism', 
                   '[H] Coenzyme transport and metabolism', 
                   '[I] Lipid transport and metabolism', 
                   '[P] Inorganic ion transport and metabolism', 
                   '[Q] Secondary metabolites biosynthesis, transport and catabolism',
                   '[X] Mobilome: prophages, transposons',
                   '[R] General function prediction only', 
                   '[S] Function unknown', 
                   '[-] No results']

    # Size parameters
    size = 0.4
    fig, ax = plt.subplots(figsize=(32, 20), facecolor=facecolor)
    fig.subplots_adjust(left=0.01, right=0.5)  # You can adjust the value as needed

    from matplotlib.colors import to_rgba

    pastel_colors = ['#62d19d', '#83c7ff', '#ffc000', '#ff8383', '#db83c6', '#8c84ff']
    green, sky, yellow, pink, violet, indigo = pastel_colors

    # Create a list of colors with varying opacity
    colors_with_opacity = \
    [to_rgba(green, alpha=0.95), to_rgba(green, alpha=0.85), to_rgba(green, alpha=0.75), to_rgba(green, alpha=0.65),
     to_rgba(green, alpha=0.55), to_rgba(green, alpha=0.45), to_rgba(green, alpha=0.35), to_rgba(green, alpha=0.25), to_rgba(green, alpha=0.15),
     to_rgba(sky, alpha=0.6),
     to_rgba(yellow, alpha=0.85), to_rgba(yellow, alpha=0.75), to_rgba(yellow, alpha=0.65), to_rgba(yellow, alpha=0.55), to_rgba(yellow, alpha=0.45),
     to_rgba(pink, alpha=0.95), to_rgba(pink, alpha=0.85), to_rgba(pink, alpha=0.75), to_rgba(pink, alpha=0.65),
     to_rgba(pink, alpha=0.55), to_rgba(pink, alpha=0.45), to_rgba(pink, alpha=0.35), to_rgba(pink, alpha=0.25),
     to_rgba(violet, alpha=0.6),
     to_rgba(indigo, alpha=0.8), to_rgba(indigo, alpha=0.6), to_rgba(indigo, alpha=0.4)]
                        
    inner_colors = pastel_colors
    outer_colors = colors_with_opacity

    if num_labels:
        ax.pie(vals,
            radius=1.4,
            startangle=140,
            colors=outer_colors,
            labels=cogs_label_nonzero_num,
            labeldistance=0.75,
            rotatelabels=True,
            textprops={'color':font_color, 'fontsize': 20, 'fontweight': 'bold'},
            wedgeprops=dict(width=size, edgecolor='w'))

    else:
        ax.pie(vals,
            radius=1.4,
            startangle=140,
            colors=outer_colors,
            labels=cogs_label_nonzero,
            labeldistance=0.85,
            rotatelabels=True,
            textprops={'color':font_color, 'fontsize': 22, 'fontweight': 'bold'},
            wedgeprops=dict(width=size, edgecolor='w'))


        ax.add_artist(con)

    ax.pie(groupsum, 
        radius=1.4-size,
        startangle=140,
        colors=inner_colors,
        # labels=groups_label,
        labeldistance=0.2,
        rotatelabels=True,
        textprops={'color':font_color, 'fontsize': 25},
        wedgeprops=dict(edgecolor='w'))
    
    ax.axis('equal')

    # Create a pie chart
    plot_title = group_name + ' ' + types + ' ' + 'Genes'
    save_title = group_name + '_' +  types + '_' +  'Genes'
    figure_name = save_title + '.png'

    # Set a title
    # ax.set_title(plot_title, fontsize=60, pad=30, color=font_color)

    # Add a legend
    fig.legend(bbox_to_anchor=(0.75, 0.5), loc='center', labels=cogs_legend, fontsize=30)
    plt.savefig(save_title + '.png', facecolor=facecolor)


if __name__ == "__main__":
    main()
