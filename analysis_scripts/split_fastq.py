from Bio import SeqIO
from pathlib import Path
import re
import sys

# Parse command-line arguments
if len(sys.argv) != 4:
    print("Usage: python3 split_fastq.py <input.fastq> <barcode_list.txt> <filecode>")
    sys.exit(1)

fastq_path = Path(sys.argv[1])
barcode_list_file = Path(sys.argv[2])
filecode = sys.argv[3]

# Read barcode list from file
barcode_list = {}
with open(barcode_list_file) as f:
    for line in f:
        sample_name, barcode = line.strip().split()
        barcode_list[barcode] = sample_name


# Create output files for each barcode
output_files = {}
for barcode in barcode_list:
    sample_name = barcode_list[barcode]
    output_files[barcode] = open(f"{sample_name}_{barcode}_{filecode}.fastq", "w")

# Scan FASTQ file for each barcode
with open(fastq_path) as handle:
    for record in SeqIO.parse(handle, "fastq"):
        s = record.description
        barcode = re.split('[:+]', s)[-2]
        if barcode in barcode_list:
            sample_name = barcode_list[barcode]
            SeqIO.write(record, output_files[barcode], "fastq")

# Close output files
for barcode in barcode_list:
    output_files[barcode].close()
