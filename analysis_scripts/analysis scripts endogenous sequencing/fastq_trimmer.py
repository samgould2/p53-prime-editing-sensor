import gzip
from string import Template
import sys
from pathlib import Path

if len(sys.argv) != 3:
    print("Usage: python3 fastq_trimmer.py <R1_FILE> <R2_FILE>")
    sys.exit(1)

R1_FILE= Path(sys.argv[1])
R2_FILE= Path(sys.argv[2])

#---Output file naming----

R1 = sys.argv[1]
r1_new = R1.split('/')[-2] + '_' + R1.split('/')[-1]
r1_new = './trimmed_fastq_files_100nt/' + r1_new.split('.')[0] + '_trimmed100.fastq'

R2 = sys.argv[2]
r2_new = R2.split('/')[-2] + '_' +R2.split('/')[-1]
r2_new = './trimmed_fastq_files_100nt/' + r2_new.split('.')[0] + '_trimmed100.fastq'


#---GLOBAL PARAMETERS-----
GZ=False
Trim_length = 100

def fastq_reader(fname, gz=False):
    _open = gzip.open if gz else open
    with _open(fname, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                identifier = line
            elif i % 4 == 1:
                read = line
            elif i % 4 == 3:
                quality = line
                yield identifier, read, quality #protospacer (index2), read, quality score


fastq_template = Template(f'$identifier\n$sequence\n$plus_sign\n$quality_score\n')
plus_sign = '+'

with open(file=r1_new, mode='a+') as file1, open(file=r2_new, mode='a+') as file2:

    for i, ((id1, r1, q1), (id2, r2, q2)) in enumerate(zip(fastq_reader(R1_FILE, gz=GZ), fastq_reader(R2_FILE, gz=GZ)), 1):

        #[:-1] removes the "\n" from the line
        id1 = id1[:-1]
        r1 = r1[:-1]
        q1 = q1[:-1]

        id2 = id2[:-1]
        r2 = r2[:-1]
        q2 = q2[:-1]
        
    
        file1.write(fastq_template.substitute(identifier=id1,
                                            sequence=r1[:Trim_length],
                                            plus_sign=plus_sign,
                                            quality_score=q1[:Trim_length]))
        
        file2.write(fastq_template.substitute(identifier=id2,
                                            sequence=r2[:Trim_length],
                                            plus_sign=plus_sign,
                                            quality_score=q2[:Trim_length]))