#!/usr/bin/env python
# prepare_reference_transcript_gtf

import numpy as np
import pandas as pd
import sys
import os
import subprocess 
import shutil
import logging 
from gtfparse import read_gtf
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) #gtfparse future warnings

logger = logging.getLogger()
srcPath = os.path.dirname(os.path.realpath(__file__)) 
sys.path.insert(1, srcPath)

'''
create_output_dir(args)
determine_order(gencode)
parse_gencode_reference(gencode_gtf, gene)
replace_geneid(df, classf)
subset_reference
parse_transcriptome_gtf(input_gtf, gene, classf)
parse_transcript(gencode, parsed_transcriptome, transcript, wobble, order)
filter_parsed_transcript(gene, gencode, parsed, output_log)
class_by_transcript_pd(transcript, All_FilteredParsed)
class_by_transcript(transcript, All_FilteredParsed)
generate_split_table(list2table,colname)
reidentify_isoform_dataset(input_file)
'''

'''
create output directories <TargetGenesRef> and <TargetGenes>
'''
def create_output_dir(args):
    # create output directories 
    args.ref_output_dir = args.output_dir + "/TargetGenesRef"
    args.gene_root_dir = args.output_dir + "/TargetGenes/" + args.genename
    args.gene_stats_dir = args.gene_root_dir + "/Stats/"                          
    args.gene_tracks_dir = args.gene_root_dir + "/Tracks/"  
                              
     
    if not os.path.exists(args.ref_output_dir):
        os.makedirs(args.ref_output_dir)
        os.makedirs(args.output_dir + "/TargetGenes/" )
            
    if os.path.isdir(args.gene_root_dir):  
        shutil.rmtree(args.gene_root_dir) 
    
    for r in [args.gene_root_dir, args.gene_stats_dir, args.gene_tracks_dir]:
        os.mkdir(r)                    
    
    return(args)


'''
**** Scenario 
sense ==> start|--exon1--|end--->---start|--exon2--|end
start, end, exon_number, transcript_id, start_end
10, 20, 1, X, 1020
40, 53, 2, X, 4053

antisense ==> start|--exon2--|end---<---start|--exon1--|end
start, end, exon_number, transcript_id, start_end
60, 50, 1, Y, 6050
40, 32, 2, Y, 4032

*** Application; use the start and end of one exon to determine order
sense: 40 - 10 > 0 
antisense: 40 - 60 < 0

By ascending start, descending end:
A |---exon 1----|end ---->---start|---exon 2----|  
B |---exon 1----|end ---->---start|---exon 2----|
C |-------------------exon 1--------------------|
read:
1. C exon 1
2. A exon 1, B exon 1
3. A exon 2, B exon 2
therefore avoids downstream complications with parse_gencode_reference given matching ends
A exon 1 and exon 2 would be assigned as different exons, despite exon 2 having a match end coordinate with C exon 1

'''
def determine_order(gencode, gene):

    # take the first transcript as example
    dat = gencode.loc[gencode["transcript_id"] == gencode["transcript_id"][0]]
    
    if dat["start"][1] - dat["start"][0] > 0:
        order = "sense"
        gencode = gencode.sort_values(['start','end'],ascending=[True, False])
    else:
        order = "antisense"
        # switch start and end coordinates of exon in antisense fashion
        gencode = gencode.rename({'start': 'end', 'end': 'start'}, axis=1) 
        # sort in descending order with exon 1 in antisense fashion
        gencode = gencode.sort_values(['start','end'],ascending=[False, True])

    gencode = gencode.reset_index(drop=True)

    return gencode, order 


    
'''
subset_reference(args)
grep genename from reference annotation using linux command
'''
def subset_reference(args):   
    print("Subsetting", args.genename, "from", args.reference)
    grepcmd = "grep {g} {f} > {d}/{o}".format(g=args.genename,f=args.reference,d=args.ref_output_dir, o=args.genename + '_gencode.gtf')
    subprocess.run(grepcmd, shell=True)
    
'''
Aim: Parse Gencode reference gtf of target gene and filter on exon genes
*Important* Check this for every target gene that the final coordinates and exon number make sense

1/ Obtain all coordinates of known exons
2/ Determine direction/order (sense or antisense)
3/ Flatten gene structure by relabelling exon number based on the start and end coordinates of exons

Scenario 1, 2, 3) Denoted the same exon with same start and/or end coordinates. 
Scenario 4 also denoted the same exon, if within the start and end coordinates are within the reference exon 
ref: ----exon------
1:   ----exon------
2:   ----exon--
3:     --exon------
4:     --exon--
The reference exon refers to the exon before   

#### Note 2:
==> antisense
situation 1:                                 outcome:
A |end2___start2|-----<------|end1___start1| (exon 2, exon 1)
B |end1______________________________start1| (exon 1)
   
situation 2:                                 outcome:
A |end2___start2|-----<------|end1___start1| (exon2, exon 1)
B |end1__________start1|                     (exon2)

using current approach with dictionary:
situation 1, all three exons are identified the same due to matching start and end coordinates 
therefore to differentiate the two exons in A:
if A_end2 = B_end1, and from dictionary A_end1 =/= B_end1 and A_start2 < A_end1, then update exon number (+1) 
=/= indirectly equate 

==> sense
situation 1:                                 outcome:
A |start1___end1|-----------|start2___end2|  (exon 1, exon 2)
B |start1_____________________________end2|  (exon 1)
to differentiate the two exons in A:
if A_end2 = B_end2 and A_start2 > A_end1, then update exon number (+1) 

   
situation 2:                                 outcome:
A |start1___end1|-----------|start2___end2|  (exon1, exon 2)
                   B |start1__________end2|  (exon2)
*prev_start = A_start1

situation 2:                                 outcome:
A |start1___end1|-----------|start2___end2|  (exon1, exon 2)
B |start1__________end1|                     (exon1)

'''
def parse_gencode_reference(args):
      
    # subsetted from gencode annotation
    gencode_gtf = args.ref_output_dir + "/" + args.genename + "_gencode.gtf"
  
    ## ----- Read and obtain coordinates
    print("**** Extracting gtf", file=sys.stdout)
    logger.disabled = True
    gencode = read_gtf(gencode_gtf)
    logger.disabled = False
    
    # obtain coordinates of known exons, note exons without exon_number in gencode are dropped
    gencode = gencode.loc[(gencode['feature'] == 'exon') & (gencode['gene_type'] == 'protein_coding') & 
                          (gencode["gene_name"] == args.genename), ["seqname","start","end","exon_number","transcript_id"]]
    gencode["start_end"] = gencode['start'].astype(str) + gencode['end'].astype(str)
    gencode.replace("", float("NaN"), inplace=True)
    
    # only using transcripts that are protein-coding (removing lncRNA, antisense to gene etc)
    gencode.dropna(subset = ["exon_number"], inplace=True)
 
    # reindex after removing empty exon_number column
    gencode.index = np.arange(0, len(gencode))
    
    ## ----- Determine direction
    # sort by start coordinates in order depending of if sense or antisense 
    gencode, order = determine_order(gencode, args.genename) 
    
    # numeric 
    gencode[["exon_number"]] = gencode[["exon_number"]].apply(pd.to_numeric)
    
    ## ----- Record alternative promoters
    # exon 1 of other transcripts
    altPromoters = gencode.loc[(gencode['exon_number'] == "1")]["start_end"].values
    
    ## ----- Record alternative terminators
    # maximum (last exon) of known transcripts
    max_exons = gencode.groupby('transcript_id')['exon_number'].max()
    altTerminators = []
    for transcript_id, max_exon in max_exons.items():
        start_end = gencode.loc[(gencode['transcript_id'] == transcript_id) & (gencode['exon_number'] == max_exon), 'start_end'].values[0]
        #print(f'Transcript ID: {transcript_id}, Max Exon: {max_exon}, Start_End: {start_end}')
        altTerminators.append(start_end)
    
    # drop duplicated exons with same start and exons to simplify the relabel
    gencode = gencode.drop_duplicates(subset=["start_end"], keep='first')
    gencode = gencode.reset_index(drop=True)
    
    # sort the exons into order from start coordinates 
    # note: to avoid issues downstream whereby one exon could be annotated to different exon number of known transcripts 
    # essentially to flatten out all known transcripts into one
    # relabel the exon_number so the exon_number goes in assending order if the current value is < previous value in the column
    prev_final = gencode.loc[0,'exon_number'] 
    prev_start = int(gencode.loc[0,'start']) 
    prev_end = int(gencode.loc[0,'end']) 
    start_dict = {}
    end_dict = {}
    list_final = []
    for index, row in gencode.iterrows():
        exon = row["exon_number"]
        start = int(row["start"])
        end = int(row["end"])

        # if the current value is the same as the previous value
        # if the start or the end coordinate is the same as before (but not situation 2) , keep the previous value 
        # otherwise update
        
        # Check various conditions and update dictionaries accordingly
        found_start = next((v for k, v in start_dict.items() if k == start), None)
        found_end = next((v for k, v in end_dict.items() if k == end), None)
        
        # matching ends - other end coordinates from identical assigned exons
        matching_ends = [key for key, value in end_dict.items() if value == found_end]

        if found_start is not None:
            list_final.append(found_start)
            start_dict[start] = found_start
            end_dict[end] = found_start
        # see NOTE 2, and note the "not any" i.e = FALSE
        elif found_end is not None and order == "antisense" and not any(start < x for x in matching_ends):
            list_final.append(found_end)
            start_dict[start] = found_end
            end_dict[end] = found_end
        elif found_end is not None and order == "sense" and not any(start > x for x in matching_ends):
            list_final.append(found_end)
            start_dict[start] = found_end
            end_dict[end] = found_end
        elif start == prev_start:
            list_final.append(prev_final)
            start_dict[start] = prev_final
            end_dict[end] = prev_final
        elif end == prev_end:
            list_final.append(prev_final)
            start_dict[start] = prev_final
            end_dict[end] = prev_final
        else:
            prev_final = int(prev_final) + 1
            list_final.append(prev_final)
            start_dict[start] = prev_final
            end_dict[end] = prev_final

        prev_start = start
        prev_end = end

    # given creating an empty list, the first entry can be exon 1 or exon 2 
    # exon 1 if the top two entries belong to exon 1, or exon 2 if the second entry belongs to exon 2
    # therefore adjust the column in the later case by minus each element by 1, so starting from exon 1
    if int(list_final[0]) != 1:
        list_final[:] = [int(num) - 1 for num in list_final]
        
    def assert_ascending_order(numbers):
        numbers = list(set(numbers))
        assert numbers, "List is empty."
        assert numbers[0] == 1, "List does not start with 1."
        assert all(a < b for a, b in zip(numbers, numbers[1:])), "List is not in ascending order starting from 1."
    
    list_final = [int(x) for x in list_final]
    assert_ascending_order(list_final)
    gencode["updated_exon_number"] = list_final
    gencode[["updated_exon_number"]] = gencode[["updated_exon_number"]].apply(pd.to_numeric)
    
    gencode['altPromoters'] = np.where(gencode['start_end'].isin(altPromoters), True, False)
    gencode["altTerminators"] = np.where(gencode['start_end'].isin(altTerminators),True,False)
    
    # write output
    gencode.to_csv(args.gene_stats_dir + args.genename + "_flattened_gencode.csv",index=False)

    return(gencode, order)


def replace_geneid(args, df):
    if df['gene_id'].str.startswith('PB.').any():
        print('Converting gtf file from SQANTI: PB gene id to gene name')
        
        # read classification file 
        print("Using", args.input_class)
        classf = pd.read_csv(args.input_class, sep = "\t")
        
        # Create a dictionary mapping PB gene IDs to associated gene names
        gene_id_name = dict(zip(["PB." + i[1] for i in classf['isoform'].str.split(".")], classf['associated_gene']))
        
        # Replace PB gene IDs in the gtf file with the corresponding dictionary values
        df['gene_id'] = df['gene_id'].replace(gene_id_name)
        
        # filter isoforms that are classified as associated to gene of interest using classfile 
        # don't just use the PBID.X from GTF as SQANTI3 misclassifies genes with the same PBID.X even though different
        df = df[df["transcript_id"].isin(list(classf[classf["associated_gene"]==args.genename]["isoform"]))]
           
    return df
    

def parse_transcriptome_gtf(args, order):
    
    '''
    Aim: Parse Input Gtf of target gene 
    '''
        
    # read gtf and filter only on exon (should only essentially remove "transcript" in gtf)
    logger.disabled = True
    df = read_gtf(args.input_gtf)
    logger.disabled = False
    
    # check if needs replacing the gene id column in the gtf from "PB.XX" to the gene name 
    df = replace_geneid(args, df)
    
    print("Total number of isoforms:", len(df["transcript_id"].unique()))
    print("**** Extracting for transcripts associated with:", args.genename)
    df = df[(df["gene_id"] == args.genename) & (df["feature"] == "exon")]
    
    if order == "antisense": 
        # reverse the rows
        df = df.iloc[::-1] 
        # rename start and end; reverse strand, the larger coordinate is the start
        df = df.rename({'start': 'end', 'end': 'start'}, axis=1)  
        
    print("Number of detected transcripts:", len(df["transcript_id"].unique()))
   
    return(df)

   
'''
Aim: parse through each transcript to classify the splicing events based on coordinates
: gencode = parsed gencode using the updated_exon_number column 
: transcript = the transcript ID 
: wobble = the maximum number of bp from the splice junction (SJ) to still be called as the same 

Output list: <transcript_id>; <transcript_exon>; <classification> <gencode_exon>

1/ Filter through the transcriptome for the relevant transcript to extract exon coordinates, re-index to 1
2/ Iterate through each exon of that transcript, capturing the start and end coordinates 
3/ Iterate through each exon of the reference against that of the target transcript, comparing the start and end coordinates

*** Exon Classification ***
** Match **
# Target exon coordinates exactly the same as the Gencode exon coordinates ==> "Match"
# Target exon coordinates within wobble range of the Gencode exon coordinates ==> "Match"

** Truncated **
# Target start coordinates exact but Target end coordinates less than wobble range ==> "TruncatedA3"
# Target end coordinates exact but Target start coordinates less than wobble range ==> "TruncatedA5"
# Target start and exon coordinates less than the wobble range ==> "TruncatedBothA3A5"

** Extended **
# Target start coordintes exact, but end coordinates more than wobble range ==> "ExtendedA3"
# Target end coordinates exact, but start coordinates more than wobble range ==> "ExtendedA5"

** Intron Retention **
# Start end coordinates exact, but end coordinates match other downstream gencode end ==> IRMatch

** None **
# No match of the start or end coordinates to the reference gencode exon
'''
def parse_transcript(gencode, parsed_transcriptome, transcript, wobble, order):

    # info to extract from the gencode reference 
    max_exon = max([int(i) for i in gencode['updated_exon_number']])
    gencode['updated_exon_number'] = gencode['updated_exon_number'].astype(str)
    list_gencode_ends = gencode["end"]
    
    # Filter 
    parsed = []
    filtered = parsed_transcriptome[parsed_transcriptome['transcript_id'] == transcript]
    
    # Reindex
    filtered.index = np.arange(1, len(filtered) + 1)
    
    # Iterature through each exon
    for index, row in filtered.iterrows():
        
        # capturing transcript start and end coordinates
        tranID = row["transcript_id"]
        tranStart = int(row["start"])
        tranEnd = int(row["end"])
        tranExon = str(index)
        
        # Match each exon coordinate of detected target transcripts to each exon in reference gencode
        for gencode_index, gencode_row in gencode.iterrows():
            
            # capturing gencode reference start and end coordinates 
            genStart = int(gencode_row["start"])
            genEnd = int(gencode_row["end"])
            genExon = str(gencode["updated_exon_number"][gencode_index])
            
            diff = str(abs(genStart-tranStart) + abs(genEnd-tranEnd))

            # For Intron retention classification 
            # Find the gencode end sites that have matching start sites to the target transcript start
            paired_ends = list(gencode[gencode["start"] == tranStart]["end"])
            # Identify the remaining gencode end sites for further intron retention classification
            other_gencode_ends = list(set(list_gencode_ends) - set(paired_ends))
            if tranEnd in other_gencode_ends: other_gencode_ends.remove(tranEnd)
            
            appText = str(tranID + ";Transcript_Exon"+ tranExon)
            noText = tranID + ";Transcript_Exon"+ tranExon +";No_Gencode_" + gencode["updated_exon_number"][gencode_index] + ";0"
            # Exact Coordinates ==> MATCH    
            if(tranStart == genStart and tranEnd == genEnd):
                parsed.append(appText +";Match_Gencode_" + genExon + ";0")
                
            # Exon Coordinates within wobble range ==> MATCH
            elif (genEnd - wobble <= tranEnd <= genEnd + wobble and 
                  genStart - wobble <= tranStart <= genStart + wobble):
                parsed.append(appText +";Match_Gencode_" + genExon + ";" + diff)
            
            # Start site the same but the end coordinates match downstream exon coordinates rather than the corresponding gencode end
            # wobble within other gencode endd
            elif (genStart - wobble <= tranStart <= genStart + wobble and
                  any(tranEnd in range(i-wobble, i+wobble) for i in other_gencode_ends)):
                IR_exon = [i for i in other_gencode_ends if tranEnd in np.arange(i-wobble,i+wobble,1)][0]
                IR_exon_df = gencode[gencode["end"] == IR_exon]
                # diff = 0 as matching with downstream exon 
                parsed.append(appText + ";IRMatch_Gencode_" + genExon + ";0")
                parsed.append(appText + ";IRMatch_Gencode_" + max(IR_exon_df["updated_exon_number"]) + ";0")

            elif order == "sense":
                '''
                        5               10
                        A5------>-------A3
                
                    A5------>-----------A3              Extended A5
                        A5------>-----------A3          Extended A3
                              A5-------->---A3          Truncated A5
                        A5---<----A3                    Truncated A3
                        A5---<----A3                    TruncatedBothA3A5
                    A5------>---------------A3          ExtendedBothA3A5       
                '''
                # Start coordinates the same, but transcript end shorter than gencode end ==> TruncatdA3
                if (genStart - wobble <= tranStart <= genStart + wobble):
                    if (tranEnd <= genEnd - wobble):
                        parsed.append(appText +";TruncatedA3_Gencode_" + genExon + ";" + diff)
                    elif (tranEnd > genEnd + wobble):
                        parsed.append(appText +";ExtendedA3_Gencode_" + genExon + ";" + diff)
                    else:
                        parsed.append(noText)
                
                elif (genEnd - wobble <= tranEnd <= genEnd + wobble):
                    if (tranStart >= genStart + wobble):
                        parsed.append(appText +";TruncatedA5_Gencode_" + genExon + ";" + diff)
                    elif (tranStart < genStart - wobble): 
                        parsed.append(appText +";ExtendedA5_Gencode_" + genExon + ";" + diff)
                    else:
                        parsed.append(noText)
                
                elif (genStart < tranEnd < genEnd - wobble and genEnd > tranStart > genStart - wobble):
                    parsed.append(appText +";TruncatedBothA3A5_Gencode_" + genExon + ";" + diff)

                elif (genEnd + wobble < tranEnd and tranStart < genStart - wobble):
                    parsed.append(appText +";ExtendedBothA3A5_Gencode_" + genExon + ";" + diff)             
                                   
                else: # No classification to the reference gencode exon
                    parsed.append(noText)


            elif order == "antisense":
                '''
                        5               10
                        A3------<-------A5
                
                        A3------<-----------A5          Extended A5
                    A3----------<-------A5              Extended A3
                        A3------<---A5                  Truncated A5
                           A3---<-------A5              Truncated A3
                    A3---<------------------A5          ExtendedBothA3A5
                           A3---<---A5                  TruncatedBothA3A5
                '''
                # A3 same
                if (genEnd - wobble <= tranEnd <= genEnd + wobble):
                    if (tranStart >= genStart + wobble):
                        parsed.append(appText +";ExtendedA5_Gencode_" + genExon + ";" + diff)                     
                    elif (tranStart <= genStart - wobble):    
                        parsed.append(appText +";TruncatedA5_Gencode_" + genExon + ";" + diff)                                         
                    else:
                         parsed.append(noText)
                    
                elif (genStart - wobble <= tranStart <= genStart + wobble):
                    if (tranEnd <= genEnd - wobble):
                        parsed.append(appText +";ExtendedA3_Gencode_" + genExon + ";" + diff)
                    elif (tranEnd >= genEnd + wobble):
                        parsed.append(appText +";TruncatedA3_Gencode_" + genExon + ";" + diff)
                    else:
                         parsed.append(noText)
                
                elif (genStart > tranEnd > genEnd + wobble and genEnd < tranStart > genStart - wobble):
                    parsed.append(appText +";TruncatedBothA3A5_Gencode_" + genExon + ";" + diff)
                    
                elif (genEnd - wobble > tranEnd > genEnd and genStart > tranStart > genStart + wobble):
                    parsed.append(appText +";ExtendedBothA3A5_Gencode_" + genExon + ";" + diff)
                    
                else: # No classification to the reference gencode exon
                    parsed.append(noText)
    
    # add empty ("0") column for classifications other than IR 
    # IR_Gencode entries have additional identifier of the direction of IR
    finalised = []
    for element in parsed:
        split_data = element.split(';')
        if len(split_data) < 5:
            split_data.append('0')
        finalised.append(';'.join(split_data))
              
    return(finalised)


 
'''
Aim: Filter the parsed output to prioritise the output 
Note: due to multiple reference transcripts with same or similar exons, 
the parse function generates redundant classifications
i.e one target exon can be complete match to one reference exon but truncated/extended to another reference exon
:gencode: parsed output gencode gtf from parse_gencode_reference()
:parsed: output from parse_transcript()


On the parsed output i.e. for each transcript, prioritise ranking for
- each gencode exon 
- each transcript exon

*** Rank Priorities ***
1. Exact match for that gencode exon
2. IR Match; start coordinates most important and so overrides Extension & Truncation
3. Extension > Truncation (like with exon 1), as longer transcript/exon is representative
4. Extension/Truncation > TruncatedBothA3A5, as at least one splice site match 
    
Note: each gencode exon stored in the reference is in the parsed file, even if no classifications recorded for it
therefore no need to know the number of gencode exons beforehand

'''
def filter_parsed_transcript(gencode, parsed, output_log):
  
    def filterRank(df,col):
        retain = []
        rank = ["Match","IRMatch","ExtendedA3","ExtendedA5","TruncatedA3","TruncatedA5","TruncatedBothA3A5","No"]
        
        for i in df[col].unique():
            subset = df[(df[col]==i)]
            if subset.shape[0] > 1:
                for r in rank:
                    if r in subset["Class"].values:
                        freq = list(subset["Class"].values).count(r)
                        
                        if freq > 1 and r != "IRMatch":
                            subsetR = subset.loc[subset["Class"] == r,]
                            mdiff = min(subsetR["Diff"].values)
                            
                            if r not in ["No"]:
                                retainElement = list(subset.loc[(subset["Class"] == r) & (subset["Diff"] == mdiff),"Original"])
                                assert len(retainElement) == 1, "Multiple gencode exons assigned: " + subset["Transcript"].values[0]
                                retain.extend(retainElement)
                                
                        else:
                            retain.extend(list(subset.loc[subset["Class"] == r,"Original"]))
                        
                        break
            else:
                retain.extend(list(df.loc[df[col]==i,"Original"]))

        retainDf = df.loc[df["Original"].isin(retain),]
        return(retainDf)
    
    parsedDf = class_by_transcript_pd(parsed[1].split(";",3)[0],parsed).drop_duplicates()
    # apply ranking to gencode exons
    GencodeParsedDf = filterRank(parsedDf,'GencodeExon')
    # apply ranking to transcript exons
    ExonParsed = filterRank(GencodeParsedDf,'TranscriptExon')
    
    # remove duplicated IRMatch entries that are not parsed
    ExonParsed = ExonParsed.drop_duplicates()
    
    assert len(ExonParsed[~ExonParsed['Class'].str.contains('IRMatch')]['TranscriptExon'].unique()) == len(ExonParsed[~ExonParsed['Class'].str.contains('IRMatch')]['TranscriptExon']), "Duplicated transcript exon entries"
    assert len(ExonParsed[~ExonParsed['Class'].str.contains('IRMatch')]['GencodeExon'].unique()) == len(ExonParsed[~ExonParsed['Class'].str.contains('IRMatch')]['GencodeExon']), "Duplicated gencode exon entries"
    
    for element in ExonParsed["Original"]:
        output_log.write(element + "\n")
    
    return(ExonParsed["Original"])


'''
Aim: Quickly filter output from filter_parsed_transcript() using the transcript ID
# Used in multiple functions downstream
:transcript = transcript ID

Output: 
:class_transcript_exon = list of transcript exons classified 
:class_gencode_exon = list of corresponding gencode exons classified with the classification types
'''
def class_by_transcript(transcript, All_FilteredParsed):
    
    # AFP = All_FilteredParsed Tuple
    # Create a tuple for filtering from the big list using the first element as transcript id and second element as class
    AFP = [(i.split(";",3)[0],i.split(";",3)[1],i.split(";",3)[2]) for i in All_FilteredParsed]
    class_transcript_exon = []
    class_gencode_exon = []
    
    for i in AFP:
        # filter the tuple based on the trnascript_id and capture the corresponding class
        if i[0] == transcript:
            class_transcript_exon.append(str(i[1]))
            class_gencode_exon.append(str(i[2]))

    return class_transcript_exon, class_gencode_exon
    

def class_by_transcript_pd(transcript, All_FilteredParsed):    
    parsed_df = pd.DataFrame([[i, 
                                 i.split(";",4)[0],
                                 i.split(";",4)[1].split("_",2)[1], 
                                 i.split(";",4)[2].split("_",2)[0], 
                                 i.split(";",4)[2].split("_",0)[0],
                                 i.split(";",4)[2].split("_",1)[1],
                                 i.split(";",4)[1],
                                 int(i.split(";",4)[3]),
                                 i.split(";",4)[4]] for i in All_FilteredParsed], 
                 columns=["Original","Transcript","TranscriptExon", "Class","ClassGencodeExon","GencodeExon","Labels","Diff","IRDir"])
    
    output_df = parsed_df.loc[parsed_df["Transcript"] == transcript,]
    output_df["GencodeExon"] = [int(i.replace("Gencode_","")) for i in output_df["GencodeExon"]]
    return(output_df)


def generate_split_table(list2table,colname):
    
    # Create as dataframe, split the column by comma into two 
    df = pd.DataFrame({'Cat':list2table})
    df[['transcript_id', colname]] = df['Cat'].str.split(',', 1, expand=True)
    df = df.drop('Cat', axis=1)
    df = df.rename(columns={"transcript_id":"transcriptID"})
    
    return(df)


'''
#### Function only applicable for merged datasets ###
Aim: Define from the isoforms of a merged dataset (i.e. PacBio and ONT) from the respective dataset
params:
    input_file = dataframe contains an "isoform" column with the PacBio and ONT isoforms merged 
    ** Important **: PacBio ID always before ONT ID if isoform present in both datasets, separated by "_"

Output: Appends two additional columns (dataset, and ONT isoforms)
'''
def reidentify_isoform_dataset(input_file):

    dataset = []
    ont_isoforms = []
    isoseq_isoforms = []
    # loop through each row of the input file
    for index, row in input_file.iterrows():
        # isoform present in both dataset with the presence of "_"
        if "_" in row["isoform"]:
            dataset.append("Both")
            isoseq_isoforms.append(row["isoform"].split("_")[0])
            ont_isoforms.append(row["isoform"].split("_")[1])
        # Iso-Seq only
        elif "PB" in row["isoform"]:
            dataset.append("Iso-Seq")
            isoseq_isoforms.append(row["isoform"])
            ont_isoforms.append("NA")
        elif "TALON" or "ENS":
            dataset.append("ONT")
            isoseq_isoforms.append("NA")
            ont_isoforms.append(row["isoform"])
        else:
            exit()

    input_file["Dataset"] = dataset 
    input_file["ONT_isoforms"] = ont_isoforms
    input_file["IsoSeq_isoforms"] = isoseq_isoforms
    
    return(input_file)
