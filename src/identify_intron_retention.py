#!/usr/bin/env python
# intron retention

import pandas as pd
import numpy as np
from prepare_and_parse import generate_split_table
from prepare_and_parse import class_by_transcript_pd

'''
Demystify the number of IR events as the number of exons with IR != IR events
:row = list of the gencode exon number classified with IR 
if row = [1, 2, 3, 5] --> Exon 1, 2, 3 with IR = 3 exons with IR but only 1 IR event 
if row = [1, 3, 4] --> same IR through Exons 3, 4 therefore total 2 IR events 

# Catalogue the number of IR events based on the series of exons
IR event is considered if the next exon is more than one up (i.e. not continuous)
'''
def demystify_IR(row):
        first = row[0] # first ordered gencode exon with IR
        prev = row[0]  # tabulate the previous exon
        count = 1 # count of IR, starting with 1
        for i in row:
            # if the exon number is not the first exon
            if(i != first):
                # if the exon number is more than the first exon + 1 i.e. a jump
                if(i != prev + 1):
                    # the include the count, and move the series along
                    count+=1
                # otherwise update the tally of previous count
                prev = i

        return(count)


'''
Aim: Tabulate the transcripts with intron retention from the All_FilteredParsed output
For the counts, apply the function to disentangle the actual number of IR events from the number of exons with IR 
1/ For transcripts with IR, order the number of exons with IR 
i.e gencode_IR = [1, 2, 3, 5] ==> Exons 1, 3 and 5 classified with IR

Last exon: if the IR is from the last exon of the updated reference

Output: 
output_df: Table of Transcripts with corresponding exon with intron retention
output_counts: The total number of transcripts with inton retention
'''        
def identify_intron_retention(args, df, All_FilteredParsed, gencode):

    max_exon = max([int(i) for i in gencode["updated_exon_number"]])

    IR = []    
    IR_Count = []
    IR_Exon1 = []
    IR_LastExon = []
    IR_ExonOverlap = []
    # Iterate through each transcript
    for transcript in df['transcript_id'].unique():
        IRMatch_Count = 0
        IROnly_Count = 0
    
        trans_df = class_by_transcript_pd(transcript, All_FilteredParsed)
        trans_df["TranscriptExon"] = [int(i.replace("Exon","")) for i in trans_df["TranscriptExon"]]
        last_exon = trans_df.loc[trans_df['TranscriptExon'].idxmax(),"GencodeExon"]
        IR_df = trans_df.loc[trans_df["Class"].isin(["IR","IRMatch"]),]
        # remove first exon IR if extend beyond A5'
        trans_df = trans_df.loc[(trans_df["GencodeExon"] != 1) & (trans_df["IRDir"] != "A5")]
        
        if len(IR_df) > 0:      
            # last exon regardless of classification 
        
            # to avoid counting Exon 1 and final exon as IR 
            #IR_df = IR_df.loc[IR_df["GencodeExon"].isin(["1",str(max_exon),str(last_exon)]) == False]
            
            if len(IR_df) > 0:
                gencode_IR = []  # list of the gencode exon number classified with intron retention
                gencode_IR = [int(item) for item in gencode_IR if item != 1]

                IR.extend(transcript + "," + row["ClassGencodeExon"] for index,row in IR_df.iterrows())
                gencode_IR.extend(int(i) for i in IR_df["GencodeExon"].values) # append the exon number

                if len(IR_df["TranscriptExon"].unique()) == 1 and IR_df["TranscriptExon"].values[0] == 1:
                    IR_Exon1.append(transcript)
                elif len(IR_df["TranscriptExon"].unique()) == 1 and IR_df["GencodeExon"].values[0] == str(max_exon):
                    IR_LastExon.append(transcript)
                elif len(IR_df["TranscriptExon"].unique()) == 1 and IR_df["GencodeExon"].values[0] == str(last_exon): 
                    IR_LastExon.append(transcript)
                else:
                    pass

                # to avoid double counting IRMatch         
                if "IRMatch" in IR_df["Class"].values: 
                    IRMatch_Count = len(np.unique(IR_df[IR_df["Class"] == "IRMatch"]["TranscriptExon"].values))
                    for t in np.unique(IR_df[IR_df["Class"] == "IRMatch"]["TranscriptExon"].values):                        
                        IRMatch_Gencode = [int(i) for i in IR_df[(IR_df["Class"] == "IRMatch") & (IR_df["TranscriptExon"] == t)]["GencodeExon"].values]
                        # + 1 to include all exons i.e exon 3 to exon 6 IR Match = 4 exons inclusive
                        IR_ExonOverlap.append(transcript + "," + str(max(IRMatch_Gencode) - min(IRMatch_Gencode)+1))
                    
                    IRMatch_GencodeAll = [int(i) for i in IR_df[IR_df["Class"] == "IRMatch"]["GencodeExon"].values]
                    gencode_IR = [i for i in gencode_IR if i not in IRMatch_GencodeAll]


                if len(gencode_IR) != 0: # if there are exons with intron retention
                    # change to integer, sort the order and apply the function
                    gencode_IR = list(set(gencode_IR))
                    gencode_IR.sort()
                    IROnly_Count = demystify_IR(gencode_IR)

                IR_Count.append(transcript + "," + str(IRMatch_Count + IROnly_Count))

    
    output_df = pd.DataFrame()
    output_counts = pd.DataFrame()
    output_exon = pd.DataFrame()
    output_overlap = pd.DataFrame()
    output_list = []
    try:
        output_df = generate_split_table(IR,"IR")
        output_df[['category', 'gencodeExon']] = output_df['IR'].str.split('_', 1, expand=True)
        output_df = output_df[['transcriptID','category','gencodeExon']]
        output_counts = generate_split_table(IR_Count,"numEvents")
        output_overlap = generate_split_table(IR_ExonOverlap,"numExonsOverlaps")
        output_list = output_df["transcriptID"].unique()

        output_df.to_csv(args.gene_stats_dir + args.genename + "_IR_transcript_level.csv", index = False)
        output_counts.to_csv(args.gene_stats_dir + args.genename + "_IR_events_counts.csv", index = False)
        output_overlap.to_csv(args.gene_stats_dir + args.genename + "_IR_numExonOverlap.csv", index = False)
    except:
        print("No transcripts with intron retention")

    
    return output_counts, output_list, IR_Exon1, IR_LastExon



