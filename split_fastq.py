import sys 
import os 
import pandas as pd 
import numpy as np
import argparse
from Bio.Seq import Seq


# Define command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="parser")  # init parser
    parser.add_argument('-i', '--fastq', help='path to input fastq file', required=True)
    parser.add_argument('-a', '--annotation', help='path to annotation file that has all primer pairs in it', required=True)
    parser.add_argument('-o', '--outdir', help='path to output directory', required=True)
    args = parser.parse_args()
    return args


# Define a function to load the primer pairs
def load_primer_pairs(args):
    primer_pairs = []  # Store all primer pairs in a list of tuples
    with open(args.annotation) as file:
        for line in file:
            splitline = line.strip().split()
            primer_id = splitline[0]
            left_primer = splitline[1]
            right_primer = splitline[2]
            primer_pairs.append((primer_id, left_primer, right_primer))
    return primer_pairs


# Define a function to check each read to see if it contains a primer pair
def find_primers(primer_pairs, seq_read):
    for primer in primer_pairs:
        primer_id = primer[0]
        primer1 = Seq(primer[1])
        primer2 = Seq(primer[2])
        primer1_rc = primer1.reverse_complement()
        primer2_rc = primer2.reverse_complement()

        primer_pair_set = [
            (str(primer1), str(primer2)),
            (str(primer1), str(primer1_rc)),
            (str(primer1), str(primer2_rc)),
            (str(primer2), str(primer1_rc)),
            (str(primer2), str(primer2_rc)),
            (str(primer1_rc), str(primer2_rc)),
        ]

        for pset in primer_pair_set:
            if (pset[0].lower() in seq_read.lower()) or (pset[1].lower() in seq_read.lower()):
                return primer_id
    return None


# Define a function to load the fastq file and read it line by line
# Define a function to load the fastq file and read it line by line
def read_fastq(args, primer_pairs):
    with open(args.fastq) as file:
        line_dict = {}
        set_id = None  # Initialize set_id (used below)
        
        for line in file:
            line = line.strip()
            
            # Check if the line starts a new record
            if line.startswith('@'):
                if len(line_dict) == 4:
                    # Process the completed record
                    seq_read = line_dict['line2']
                    set_id = find_primers(primer_pairs, seq_read)
                    
                    if set_id is not None:
                        if args.outdir.endswith('/'):
                            out_dir = args.outdir
                        else:
                            out_dir = args.outdir + "/"
                        
                        base_name = args.fastq.split("/")[-1]
                        outfile = f"{out_dir}{set_id}_{base_name}"

                        with open(outfile, 'a') as ofile:
                            lines_to_add = [line if line.endswith('\n') else line + '\n' for line in line_dict.values()]
                            ofile.writelines(lines_to_add)
                
                # Reset line_dict for the new record
                line_dict = {'line1': line}
                continue
            
            # Add subsequent lines to the dictionary
            if len(line_dict) == 1:
                line_dict['line2'] = line
            elif len(line_dict) == 2:
                line_dict['line3'] = line
            elif len(line_dict) == 3:
                line_dict['line4'] = line
            
            # If the dictionary grows too large, log an error
            if len(line_dict) > 4:
                print("line_dict is too long, is the fastq malformed?")
                break

        # Process the final record if the file ends with a complete record
        if len(line_dict) == 4:
            seq_read = line_dict['line2']
            set_id = find_primers(primer_pairs, seq_read)

            if set_id is not None:
                if args.outdir.endswith('/'):
                    out_dir = args.outdir
                else:
                    out_dir = args.outdir + "/"
                
                base_name = args.fastq.split("/")[-1]
                outfile = f"{out_dir}{set_id}_{base_name}"

                with open(outfile, 'a') as ofile:
                    lines_to_add = [line if line.endswith('\n') else line + '\n' for line in line_dict.values()]
                    ofile.writelines(lines_to_add)



# Define the main function
def main():
    args = parse_arguments()
    primer_pairs = load_primer_pairs(args)
    read_fastq(args, primer_pairs)


if __name__ == '__main__':
    main()
