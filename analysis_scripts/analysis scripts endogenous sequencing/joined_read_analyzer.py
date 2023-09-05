import gzip
import pandas as pd
import numpy as np
import sys
from pathlib import Path

if len(sys.argv) != 4:
    print("Usage: python3 joined_read_analyzer.py <joined_file> <F_HANDLE> <R_HANDLE>")
    sys.exit(1)


#importing info
joined_file = Path(sys.argv[1])
F_HANDLE = str(sys.argv[2])
R_HANDLE = str(sys.argv[3])

#global parameters
GZ=False 
MIN_QUALITY = 30

def fastq_reader(fname, gz=False):
    _open = gzip.open if gz else open
    proc_read = (lambda line: line.strip().decode()) if gz else (lambda line: line.strip())
    proc_qual = (lambda line: line.strip()) if gz else (lambda line: line.strip().encode('utf-8'))
    with _open(fname, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                identifier = proc_read(line)
            elif i % 4 == 1:
                read = proc_read(line)
            elif i % 4 == 3:
                quality = proc_qual(line)
                yield identifier, read, quality #protospacer (index2), read, quality score


high_qual_reads = []
low_quality_count = 0
unaligned = 0
for i, ((id1, r1, q1)) in enumerate(fastq_reader(joined_file, gz=GZ), 1):

    #quality score filtering
    qual = (np.frombuffer(q1, dtype=np.uint8) - 33).mean()
    
    if qual >= MIN_QUALITY:
        start = r1.find(F_HANDLE) + len(F_HANDLE)
        end = r1.find(R_HANDLE)

        if (start==-1) or (end ==-1):
            unaligned +=1

        high_qual_reads.append(r1[start:end])

    else:
        low_quality_count+=1

#generate the data frames
u, c = np.unique(high_qual_reads, return_counts=True)
df_edits = pd.DataFrame(dict(zip(['Sequence','Count'], [u,c]))).sort_values(by='Count', ascending=False).reset_index().drop(columns='index')


count1 = [len(high_qual_reads), unaligned, low_quality_count]
row_labels = ['Processed_reads', 'Unaligned', 'Low_quality']
col_labels = ['Classification', 'Count']
df_classif = pd.DataFrame(dict(zip(col_labels, [row_labels, count1])))


#and save it to the appropriate location in a dataframe after processing
s = str(sys.argv[1])
s2 = s.split('/')[-1]
jf_name = s2.split('.')[0]

file_path1 = './230801San/processed_dataframes/' + jf_name + '.csv'
df_edits.to_csv(file_path1, index=False)


file_path2 = './230801San/processed_dataframes/classification_dataframes/' + jf_name + '_read_classification.csv'
df_classif.to_csv(file_path2, index=False)