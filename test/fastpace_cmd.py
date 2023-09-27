import argparse, os, json
try:
    import fastpace
except ImportError:
    print("fastpace module is not installed. Please install it by running 'pip install fastpace'.")

def get_parsed_args():
    parser = argparse.ArgumentParser(description='Parse command arguments for a fasta file processing script.')
    parser.add_argument('--input_file', type=str, help='The path to the input fasta file', required=True)
    parser.add_argument('--output_file', type=str, help='The path to the output fasta file', required=True)
    parser.add_argument('--draw_logo', action='store_true', help='Draw sequence log of the best peptide or the peptide passed by --sequence.')
    parser.add_argument('--sequence', type=str, help='Draw sequence logo of this peptide.')
    parser.add_argument('--num_reruns', type=int, help='number of times to rerun the algorithm before returning the results.', default=1)
    parser.add_argument('--refine', type=int, help='Flag to run the refinement.', default=1)

    args = parser.parse_args()
    return args

def read_fasta(input_file):
    # Open the input file for reading
    with open(input_file, 'r') as f:
        # Read the contents of the file into a list of lines
        lines = f.readlines()

    # Initialize the header and sequence strings
    header = ""
    sequence = ""

    # Initialize a list to hold the fasta records
    fasta_records = []

    # Loop through the lines of the input file
    for line in lines:
        # If the line starts with a '>', it's a header
        if line[0] == ">":
            # If the header and sequence strings are not empty, add them to the list of fasta records
            if header and sequence:
                fasta_records.append((header.strip(), sequence.strip()))

            # Reset the header and sequence strings
            header = line
            sequence = ""
        # If it's not a header, it's a sequence line
        else:
            # Add the line to the sequence string
            sequence += line.strip()

    # Add the last header and sequence to the list of fasta records
    fasta_records.append((header.strip(), sequence.strip()))

    return fasta_records

def write_fasta(fasta_records, aligned_sequences, output_file):
    # Open the output file for writing
    with open(output_file, 'w') as f:
        # Loop through the list of fasta records
        for header, sequence in fasta_records:
            # Write the header and sequence to the output file
            f.write(header + "\n" + aligned_sequences[sequence] + "\n")

def plt_logo(scores_dict, output_path):
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import logomaker

    cols = set()
    for residue in scores_dict:
        cols = cols.union(scores_dict[residue].keys())
    cols_list = list(cols)
    
    data = {}
    for residue in scores_dict:
        row = []
        for letter in cols_list:
            if letter in scores_dict[residue]:
                row.append(scores_dict[residue][letter])
            else:
                row.append(0)
        data[int(residue)+1] = row
        
    
    logo_df = pd.DataFrame.from_dict(data, orient='index', columns=cols_list)
    
    #set styles
    sns.set(font_scale = 1.6)
    sns.set_style('ticks')
    # create Logo object
    fig, ax = plt.subplots(num='16X8',figsize=[16,8])
    logo = logomaker.Logo(logo_df,
                ax=ax,
                color_scheme='NajafabadiEtAl2017',
                vpad=.1,
                width=.8,
                shade_below=.5,
                fade_below=.5)

    logo.ax.set_ylabel('Score', fontsize=30)
    logo.ax.set_xlim([0, len(logo_df)+1])
    logo.ax.set_xticks([])
    plt.savefig(os.path.splitext(output_path)[0] +'.png', facecolor='w', edgecolor='none')
    plt.clf()


if __name__ == "__main__":
    args = get_parsed_args()
    fasta_records = read_fasta(args.input_file)
    sequences = [record[1] for record in fasta_records]
    
    if args.num_reruns > 1:
        import re
        result_json = fastpace.run_motif_discovery(sequences)
        for i in range(args.num_reruns-1):
            best_consensus = result_json['consensus']['best_motif'][2:-2]
            mask = 'X' * len(best_consensus)
            fasta_records = [(header, re.sub(best_consensus, mask, sequence)) for header, sequence in fasta_records]
            masked_sequences = [record[1] for record in fasta_records]
            result_json = fastpace.rerun_motif_discovery(sequences, masked_sequences)
    else:
        result_json = fastpace.run_motif_discovery(sequences, refine=args.refine)
    
    write_fasta(fasta_records, result_json['alignment']['aligned_sequences'], args.output_file)
    
    if args.draw_logo:
        if args.sequence:
            if args.sequence in result_json['peptides']:
                plt_logo(result_json['peptides'][args.sequence]['similarity_matrix'], args.output_file)
            else:
                print("Can't draw sequence logo for a sequence not in the dataset!")
        else:
            alignment_template = result_json['alignment']['template']
            plt_logo(result_json['peptides'][alignment_template]['similarity_matrix'], args.output_file)
    
    with open(os.path.splitext(args.output_file)[0] +'.json', 'w') as outfile:
        json.dump(result_json, outfile)
