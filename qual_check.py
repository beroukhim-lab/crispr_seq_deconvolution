import sys 
import os 
import pandas as pd 
import numpy as np
import argparse
import re
import itertools
import seaborn as sns 
import matplotlib.pyplot as plt

#Define the command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="parser")  # init parser
    parser.add_argument('-i', '--input_dir', help='path to input fastq files', required=True)
    parser.add_argument('-a', '--annotation', help='path to annotation file that has all primer pairs in it', required=True)
    args = parser.parse_args()
    return args

#Get a set of primer IDs
def get_primer_ids(args):
    with open(args.annotation) as file:
        primer_ids = []
        for line in file:
            line = line.strip().split()
            primer_id = line[0]
            primer_ids.append(primer_id)
        
        pid_set = set(primer_ids)
        return pid_set
    

#Look in input_dir and find all fastq files that begin with one of the primer IDs and end with fastq
def find_fastq(args, pid_set):
    all_files = os.listdir(args.input_dir)
    all_fastq = [fastq for fastq in all_files if fastq.endswith(".fastq")]
    all_matching_fastq = [fastq for fastq in all_fastq if any(fastq.startswith(prefix) for prefix in pid_set)]

    #Now strip the primer name from the files and return just the base fastq name
    all_fastq_basenames = []
    for fastq in all_matching_fastq:
        for pid in pid_set:
            pid_string = pid + "_"
            if fastq.startswith(pid_string):
                stripped_name = re.sub(pid_string, "", fastq)
                all_fastq_basenames.append(stripped_name)
    
    basename_set = set(all_fastq_basenames)
    return basename_set


#Now construct the read count matrix
def construct_read_count_matrix(args, pid_set, basename_set):
    pair_set = set(itertools.product(pid_set, basename_set))
    pair_string = [pair[0] + "_" + pair[1] for pair in pair_set] #Each value in pair_string is a possible file name

    #Now loop through all possible files and count the number of reads
    entries_list = []
    for pair in pair_string:
        if args.input_dir.endswith("/"):
            complete_path = args.input_dir + pair
        else:
            complete_path = args.input_dir + "/" + pair

        if os.path.exists(complete_path):
            with open(complete_path, "r") as file:
                lc = sum(1 for line in file)

                if lc%4 != 0:
                    raise ValueError(f"Malformed FASTQ file: {complete_path}. Line count not divisible by 4.")

                else:
                    num_reads = lc // 4
                    entries_list.append(int(num_reads))
        else:
            entries_list.append(0)

    
    #Now bind everything together in a single pandas df
    df_for_plotting = pd.DataFrame(pair_set)
    df_for_plotting['counts'] = entries_list
    new_names = ['primer_set', 'fastq', 'read_counts']
    df_for_plotting.columns = new_names

    #Write the read count summary matrix
    if args.input_dir.endswith("/"):
        outfile = args.input_dir + 'read_count_summary_matrix.txt'
    else:
        outfile = args.input_dir + "/read_count_summary_matrix.txt" 
    df_for_plotting.to_csv(outfile, sep = "\t", index = False)

    #Also return the output for plotting
    return df_for_plotting


#Now define a function to plot the heatmap
def plot_heatmap(args, df_for_plotting):
    plt.figure(figsize=(13,16))
    heatmap_data = df_for_plotting.pivot(index='fastq', columns='primer_set', values='read_counts')
    sns.heatmap(heatmap_data, cmap="coolwarm", cbar=True)    
    plt.title("read counts")

    #Output
    if args.input_dir.endswith("/"):
        outfile = args.input_dir + 'read_count_heatmap.pdf'
    else:
        outfile = args.input_dir + "/read_count_heatmap.pdf" 

    #Write heatmap
    plt.rcParams['pdf.fonttype'] = 42
    plt.savefig(outfile, format='pdf', bbox_inches="tight", pad_inches=0.5)



#Define the main function
def main():
    args = parse_arguments()
    pid_set = get_primer_ids(args)
    basename_set = find_fastq(args, pid_set)
    df_for_plotting = construct_read_count_matrix(args, pid_set, basename_set)
    plot_heatmap(args, df_for_plotting)



#Do the thing
if __name__ == '__main__':
    main()


