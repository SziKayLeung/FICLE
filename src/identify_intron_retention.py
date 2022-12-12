#!/usr/bin/env python
# intron retention

import pandas as pd
from .prepare_and_parse import generate_split_table
from .prepare_and_parse import class_by_transcript_pd

def identify_intron_retention(df, All_FilteredParsed, gencode):
    
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
    
    def demystify_IR(row):
        '''
        Demystify the number of IR events as the number of exons with IR != IR events
        :row = list of the gencode exon number classified with IR 
        if row = [1, 2, 3, 5] --> Exon 1, 2, 3 with IR = 3 exons with IR but only 1 IR event 
        if row = [1, 3, 4] --> same IR through Exons 3, 4 therefore total 2 IR events 
        
        # Catalogue the number of IR events based on the series of exons
        IR event is considered if the next exon is more than one up (i.e. not continuous)
        '''
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
    
    max_exon = max([int(i) for i in gencode["updated_exon_number"]])
    
    IR = []    
    IR_Count = []
    IR_Exon1 = []
    IR_LastExon = []
    IRMatch_Count = 0
    IROnly_Count = 0
    # Iterate through each transcript
    #for transcript in ["PB.8560.273_TALONT001396616"]:
    for transcript in df['transcript_id'].unique():
        trans_df = class_by_transcript_pd(transcript, All_FilteredParsed)
        
        try: 
            gencode_IR = []  # list of the gencode exon number classified with intron retention
            IR_df = trans_df.loc[trans_df["Class"].isin(["IR","IRMatch"]),]
            IR.extend(transcript + "," + row["ClassGencodeExon"] for index,row in IR_df.iterrows())
            gencode_IR.extend(int(i.split("_")[1]) for i in IR_df["GencodeExon"].values) # append the exon number

            # last exon regardless of classification 
            IR_df.loc[:,("TranscriptExon")] = [int(i.replace("Exon","")) for i in IR_df.loc[:,("TranscriptExon")]]
            last_exon = IR_df.loc[IR_df['TranscriptExon'].idxmax(),"GencodeExon"].replace("Gencode_","")
            #print(last_exon)

            if len(IR_df["TranscriptExon"].unique()) == 1 and IR_df["TranscriptExon"].values[0] == 1:
                IR_Exon1.append(transcript)
            elif len(IR_df["TranscriptExon"].unique()) == 1 and IR_df["GencodeExon"].values[0] == "Gencode_" + str(max_exon):
                IR_LastExon.append(transcript)
            elif len(IR_df["TranscriptExon"].unique()) == 1 and IR_df["GencodeExon"].values[0] == "Gencode_" + str(last_exon): 
                IR_LastExon.append(transcript)
            else:
                pass
                   
            # to avoid counting Exon 1 and final exon as IR 
            IR = [i for i in IR if i.split(",",2)[1] not in ["IR_Gencode_1","IR_Gencode_" + str(max_exon),"IR_Gencode_" + str(last_exon)]]
            
            # to avoid double counting IRMatch         
            if "IRMatch" in IR_df["Class"].values: 
                IRMatch_Count = len(np.unique(IR_df[IR_df["Class"] == "IRMatch"]["TranscriptExon"].values))
                IRMatch_Gencode = [int(i.replace("Gencode_","")) for i in IR_df[IR_df["Class"] == "IRMatch"]["GencodeExon"].values]

                gencode_IR = [i for i in gencode_IR if i not in IRMatch_Gencode]

            if len(gencode_IR) != 0: # if there are exons with intron retention
                # change to integer, sort the order and apply the function
                gencode_IR = [int(item) for item in gencode_IR if item != 1]
                gencode_IR = list(set(gencode_IR))
                gencode_IR.sort()
                IROnly_Count = demystify_IR(gencode_IR)

            IR_Count.append(transcript + "," + str(IRMatch_Count + IROnly_Count))
        except:
            pass
    
    output_df = pd.DataFrame()
    output_counts = pd.DataFrame()
    output_list = []
    try:
        output_df = generate_split_table(IR,"IR")
        output_counts = generate_split_table(IR_Count,"IR")
        output_list = output_df["transcript_id"].unique()
    except:
        print("No transcripts with intron retention")
    
    return output_df, output_counts, output_list, IR_Exon1, IR_LastExon


