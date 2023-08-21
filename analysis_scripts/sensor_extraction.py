#importing necessary packages
from collections import defaultdict
import gzip
from Bio.Align import PairwiseAligner
import numpy as np
#import matplotlib.pyplot as plt
import Bio.Seq
import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns
import Bio.Seq
import Bio.SeqIO
from Bio.SeqRecord import SeqRecord
import os
import sys
from pathlib import Path


#-------Loading in system args
# Parse command-line arguments to read in files of interest
if len(sys.argv) != 7:
    print("Usage: python3 sensor_extraction.py <input_df> <unique_protos> <R1_FILE> <R2_FILE> -o <folder_name>")
    sys.exit(1)

#example usage

input_df = pd.read_csv(Path(sys.argv[1]))
unique_protos = np.load(Path(sys.argv[2]),allow_pickle=True)

R1_FILE= Path(sys.argv[3])
R2_FILE= Path(sys.argv[4])

folder_name = str(sys.argv[6])



#----------GLOBAL VARIABLES

GZ = False #True when files are gzipped

# Minimum average read quality
MIN_QUALITY = 30 #originally set to 30

R1_CUTOFF = 15 # R1_anchor alignment score cutoff
R2_CUTOFF = 15 # R2_anchor alignment score cutoff

SPACER_CUTOFF = 16 # Spacer alignment score cutoff
#EXTENSION_CUTOFF = Needs to be dependent on the size of the 3' extension...

#R1_ANCHOR = 'CTTGTGGAAAGGACGAAACACC'  # Left sequence to anchor to
#R2_ANCHOR = 'CAACCAACGCTATCGGCATGG'  # Right sequence to anchor to (complementary strand)
R1_ANCHOR = 'TTGTGGAAAGGACGAAACACC'
R2_ANCHOR = 'CTAGCGTTCGAGTTAGGAATT'

tevopreq1_stop = 'CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAATTTTTTT' #may need to adjust the exact size of this...

#tevopreq1_stop_rev = str(Bio.Seq.Seq('CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAATTTTTTT').reverse_complement()) #may need to adjust size of this as welll...
tevopreq1_stop_rev = str(Bio.Seq.Seq('GCGTTAAACCAACTAGAATTTTTTT').reverse_complement()) #may need to adjust size of this as welll...


TEVO_CUTOFF = len(tevopreq1_stop)-6
#SCAF_CUTOFF = len(scaff_seq)-6

MIN_SPACER_LENGTH = 16 # Minimum length of spacer
MIN_SITE_LENGTH = 200 # Minimum length of edit site & tevo stop region...

# Alignment parameters (Needleman-Wunsch)
MATCH_SCORE = 1
MISMATCH_SCORE = -0.5
OPEN_GAP_SCORE = -5
EXTEND_GAP_SCORE = -0.1
TARGET_END_GAP_SCORE = 0
QUERY_END_GAP_SCORE = 0



#THIS FUNCTION MAY NEED TO BE UPDATED DEPENDING ON HOW THE DATA IS STRUCTURED
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
                yield identifier[-20:], read, quality #protospacer (index2), read, quality score

def make_aligner():
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = MATCH_SCORE
    aligner.mismatch_score = MISMATCH_SCORE
    aligner.open_gap_score = OPEN_GAP_SCORE
    aligner.extend_gap_score = EXTEND_GAP_SCORE
    aligner.target_end_gap_score = TARGET_END_GAP_SCORE
    aligner.query_end_gap_score = QUERY_END_GAP_SCORE
    return aligner

def quality_checker(q1, q2):
    """ 
    Checks quality score of the reads
    """

    qual1 = (np.frombuffer(q1, dtype=np.uint8) - 33).mean()
    qual2 = (np.frombuffer(q2, dtype=np.uint8) - 33).mean()

    low_qual = []

    if qual1 < MIN_QUALITY:
        low_qual.append('1')
    if qual2 < MIN_QUALITY:
        low_qual.append('2')
    
    if len(low_qual)==0:
        quality_out = 'good quality'
    else:
        quality_out = 'low_qual_r' + ''.join(low_qual)

    return quality_out

aligner = make_aligner()


def extension_identifier(r2, ext_sequences, tevopreq1_stop):
    """ 
    Determines ext_idx based on R2 read.
    Returns: (1) out = str describing
    (2) ext_idx = identifier of pegRNA (corresponds with pegRNA_idx, or row # with 0 start)
    """
    #r2 corresponds with extension
    #r3 = protospacer

    #initialize placeholder values
    
    out=None
    ext_idx = None

    #-----------first find the 3' extension using the tevopreq1 stop/anchor-----------

    #find the anchor
    offset_extension_end = r2.find(tevopreq1_stop)
    
    #if no perfect match, try alignment to find it (can get rid of this method for simplicity since it usually doesn't work...)
    if offset_extension_end < len(tevopreq1_stop):
        try:
            alignment_2 = sorted(aligner.align(r2, tevopreq1_stop))[0]
            score_2 = alignment_2.score
            offset_extension_end = alignment_2.aligned[0][-1][0] #last index refers to start point...
        
        except IndexError:
            out = 'unaligned tevo'
            return out, ext_idx

        if (score_2 < TEVO_CUTOFF):
            out = 'unaligned tevo'
            return out, ext_idx

        else:
            offset_extension_end = offset_extension_end
        
    else:
        offset_extension_end = offset_extension_end

    
    #---------then extract the extension sequence and identify it -------------
    extension_seq = r2[:offset_extension_end]

    #finds all of the matches to the extension
    ext_idx = np.where(ext_sequences == extension_seq)[0]

    if len(ext_idx)==0: 
        out = 'no extension match'        

    else:
        out = 'extension match'
    
    return out, ext_idx
    

def proto_extension_matcher(ext_idx, unique_protos, proto1, input_df):
    """
    Determines if the protospacer sequence matches the extension. 
    Returns:
    (1) out = str describing matching outcome
    (2) pegRNA_idx = pegRNA_idx that identifies the pegRNA of interest

    Note: due to an error, there are some duplicate pegRNAs. These should be removed, or dealt with in post-processing. 
    """
    
    pegRNA_idx=None
    #------- now identify the protospacer; this step assumes that it's the correct size (20 bp) ----------    
    proto_index = np.where(unique_protos==proto1)[0]

    if len(proto_index)==0: 
        out = 'no protospacer match'
        return out, pegRNA_idx
    
    #-------- now identify the pegRNA or identify decoupling events ---------------
    if len(ext_idx)==1:
        #should be vast majority of cases
        proto_idx_ext = input_df[input_df['pegRNA_idx'].isin(ext_idx)]['unique_proto_idx'].values[0] #this is the set of protospacers that come from the 3' extension matches

        if proto_idx_ext==proto_index:
            out = 'correct id'
            pegRNA_idx = ext_idx[0]
        else:
            out = 'decoupled extension-protospacer'

    
    elif len(ext_idx)>1:
        #in this case resort to a for loop
        proto_idx_ext = input_df[input_df['pegRNA_idx'].isin(ext_idx)]['unique_proto_idx'].values

        
        matching_set = []
        for i, val in enumerate(proto_idx_ext): #iterate through and check for matches
            if proto_index == val:
                matching_set.append(ext_idx[i])
            else:
                continue


        if len(matching_set)==0:
            out = 'decoupled extension-protospacer'

        else: 
            out = 'correct id'
            pegRNA_idx = matching_set[0] #only add count to first identified (can fix this in post-processing to account for duplicates)
        


    return out, pegRNA_idx


def sensor_extract(r1, tevopreq1_stop_rev):

    offset_sensor_end = r1.find(tevopreq1_stop_rev)

    sensor_seq = None
    
    #if no perfect match, try alignment to find it (can get rid of this method for simplicity since it usually doesn't work...)
    if offset_sensor_end < 0:
        edit_outcome = 'unaligned'

        
    else:
        #---------then extract the extension sequence and identify it -------------
        edit_outcome = 'extracted'
        sensor_seq = r1[:offset_sensor_end]


    return sensor_seq, edit_outcome


def recombination_checker(sensor_seq, sensor_handles, handle_true):

    recomb=False

    #check for recombination
    handle = sensor_seq[0:5] + sensor_seq[-5:]

    if handle==handle_true:
        recomb=False
    else:
        if handle in sensor_handles:
            recomb=True
        else:
            recomb=False
    
    return recomb

'''Prepare fastq records.'''
def to_IOSeq_rec(sensor_seq, i, q1):
    qname = 'read_' + str(i)
    record = SeqRecord(Bio.Seq.Seq(sensor_seq), id=qname, name=qname, description='', dbxrefs=[])

    #add quality score to the record
    qual = (np.frombuffer(q1, dtype=np.uint8) - 33)
    record.letter_annotations["phred_quality"] = qual[0:len(sensor_seq)]
    
    return record

def extraction_filtration(folder_name, input_df, unique_protos, R1_FILE,R2_FILE, breakpoint=False):

    """
    Takes in reads and returns dataframe containing pegRNA counts AND sensor outcomes.

    Also returns summary of outcomes including:
    - low quality (and which of the reads are low quality (or if all are low quality))
    - no extension match
    - no protospacer match
    - decoupled extension-protospacer
    - correct identification
    """

    os.mkdir(folder_name)#make the new directory (needs to be a non-existent folder)

    #write 28000 empty fastq files to a new folder
    for i in list(input_df['peg_id']):
        f_name = i + '.fastq'
        with open(folder_name + '/' + f_name, 'w') as fp:
            pass


    #also an array of the extension sequences
    ext_sequences = np.array(list(input_df['PBS_RTT_5to3']))
    
    sensor_handles = np.unique(input_df['sensor_handle'])

    
    #------initialize a dataframe for holding the pegRNA counts and sensor outcomes
    d1 = pd.DataFrame(dict(zip(['pegRNA_idx'], [list(input_df['pegRNA_idx'])])))
    cols = ['total', 'unaligned', 'extracted', 'recombination']
    z = np.zeros((len(cols), len(input_df)))
    d2 = pd.DataFrame(dict(zip(cols, z)))
    df1 = pd.concat((d1, d2), axis=1).set_index('pegRNA_idx')

    #------initialize a dataframe for holding metadata about the identification of sensors
    outcomes = ['good quality', 'low_qual_r1', 'low_qual_r2', 'low_qual_r12', 
                'no extension match', 'unaligned tevo', 'no protospacer match', 'decoupled extension-protospacer', 
                'correct id']
    outcomes_count = np.zeros(len(outcomes))
    class_df = pd.DataFrame(dict(zip(['classification', 'count'],[outcomes, outcomes_count]))).set_index('classification')

    #iterating through the reads...
    for i, ((proto1, r1, q1), (proto2, r2, q2)) in enumerate(zip(fastq_reader(R1_FILE, gz=GZ), fastq_reader(R2_FILE, gz=GZ)), 1):
        #r1 = sensor read
        #r2 = extension read
        #proto1/2 = protospacer read (index2)

        #first check the quality
        quality_out = quality_checker(q1,q2)
        class_df.loc[quality_out, 'count']+=1

        #if quality is poor, abort
        if quality_out != 'good quality':
            continue


        #else, if quality is good, continue
        else:
            out1, ext_idx = extension_identifier(r2, ext_sequences, tevopreq1_stop)
        
            if out1=='extension match':
                out, pegRNA_idx = proto_extension_matcher(ext_idx, unique_protos, proto1, input_df)
                class_df.loc[out, 'count']+=1
                #sensor classification
                #and then look at the sensor classification; probably need to put limitations on length of the read here...
                if out == 'correct id':
                    df1.loc[pegRNA_idx, 'total'] += 1
                    
                    sensor_seq, edit_outcome = sensor_extract(r1, tevopreq1_stop_rev)
                    
                    #record the edit outcome
                    if edit_outcome=='unaligned':
                        df1.loc[pegRNA_idx, edit_outcome] +=1

                    else:
                        handle_true = input_df[input_df['pegRNA_idx']==pegRNA_idx]['sensor_handle'].values[0]

                        #check for recombination (if recomb, returns true)
                        recomb = recombination_checker(sensor_seq, sensor_handles, handle_true)

                        if recomb==True:
                            df1.loc[pegRNA_idx, 'recombination'] +=1

                        else:
                            #add it to the appropriate 
                            df1.loc[pegRNA_idx, 'extracted'] +=1

                            #add the sensor read to the appropriate fastq file (numbered according to the pegRNA index)
                            out_file = folder_name + '/' + 'peg_' + str(pegRNA_idx) + '.fastq'

                            record = to_IOSeq_rec(sensor_seq, i, q1)

                            with open(out_file, 'a') as fq:
                                #and write it to the appropriate file
                                Bio.SeqIO.write(record, fq, 'fastq')

                else:
                    continue

            else:
                class_df.loc[out1, 'count']+=1
              


        if breakpoint != False:
            if i>breakpoint:
                break

    
    #and then prune the duplicates...
    #both the files and the counts table...

   
    return df1, class_df


#----------ACTUALLY RUNNING THE PIPELINE BELOW-------------

#input_df = pd.read_csv('/Volumes/sanchezrivera/samgould/p53_library_NO_DUPLICATES_with_handles.csv')
#unique_protos = np.load('unique_protos.npy',allow_pickle=True)

count_df, class_df = extraction_filtration(folder_name, input_df, unique_protos, R1_FILE,R2_FILE, breakpoint=False)

#save the counts and classification dataframes
count_df.to_csv(folder_name + '_count_df.csv')
class_df.to_csv(folder_name + '_classification_df.csv')

