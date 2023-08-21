import pandas as pd
import os
import sys
from pathlib import Path


#-------Loading in system args
# Parse command-line arguments to read in files of interest
if len(sys.argv) != 3:
    print("Usage: python3 crispresso_analysis.py <input_df> <sample_name>")
    sys.exit(1)

#example usage
p53 = pd.read_csv(Path(sys.argv[1]))
sample_name = Path(sys.argv[2])

#--------run it for all of the pegRNAs
for i in range(len(p53)):
    peg_id = p53.iloc[i]['peg_id']
    wt = p53.iloc[i]['sensor_WT_rc']
    edited = p53.iloc[i]['sensor_CorrectEdit_rc']
    proto = p53.iloc[i]['protospacer']
    wt_qwc = p53.iloc[i]['wt_qwc']
    edit_qwc = p53.iloc[i]['edit_qwc']

    #generating crispresso command
    f_path = f"./{sample_name}/{peg_id}.fastq"
    start = f"CRISPResso --fastq_r1 {f_path} "
    #removed --discard_indel_reads; leaving this in there...
    p2 = f"--amplicon_seq {wt} --expected_hdr_amplicon_seq {edited} --guide_seq {proto} --quantification_window_coordinates {wt_qwc},{edit_qwc} --suppress_report --suppress_plots --plot_window_size 15 --exclude_bp_from_left 0 --exclude_bp_from_right 0" #--discard_indel_reads

    #output folder = crispresso_output
    out_path = f" -o ./crispresso_output/{sample_name}" 

    command = start+ p2 + out_path

    #running the crispresso command
    #apparently this is bad practice (but subprocess module didn't work...)
    os.system(command)