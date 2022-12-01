#!/usr/bin/env python
# A5A3 and AF

import pandas as pd
from .prepare_and_parse import class_by_transcript
from .prepare_and_parse import class_by_transcript_pd
from .prepare_and_parse import generate_split_table

def identify_A5A3_transcripts(df, All_FilteredParsed):
    
    '''
    Aim: Tabulate the transcripts with extended or truncated exons from the All_FilteredParsed output
    i.e. Alternative 5' start and 3' end sites. 
    
    Output: 
    output_df: Table of Transcripts with corresponding exon with extension or truncation
    output_pertrans_counts: The number of extended or truncated exons per transcript
    output_counts: The total number of transcripts with extended or truncated exons
    
    '''
        
    output = []
    for transcript in df["transcript_id"].unique():
            class_transcript_exon,class_gencode_exon = class_by_transcript(transcript, All_FilteredParsed)
            # Filter through the parsed output for "Truncated" and "Extended"  
            A5A3_Exons = list(filter(lambda x: ("Truncated" in x), class_gencode_exon)) + list(filter(lambda x: ("Extended" in x), class_gencode_exon)) 
            output.append([transcript + "," + s for s in A5A3_Exons])

    # Aggregate all list to one big list
    output_df = pd.DataFrame()
    output_pertrans_counts = pd.DataFrame()
    output_counts = pd.DataFrame()
    output_transcripts = []
    try:
        output_lst = sum(output, [])
        output_df = generate_split_table(output_lst,"A5A3")[["transcript_id","A5A3"]]
        output_transcripts = list(output_df["transcript_id"].unique())

        # Output Stats
        output_df[['cate', 'gencode_exon']] = output_df['A5A3'].str.split('_', 1, expand=True)
        output_pertrans_counts = output_df.groupby(['cate','transcript_id']).size().reset_index()
        output_counts = output_df.groupby("cate").size().reset_index()
    except:
        print("No Transcripts with A5A3 truncation or extension")
    
    return output_df, output_pertrans_counts, output_counts, output_transcripts

    
def identify_alternative_first(df, All_FilteredParsed):
    
    '''
    Aim: Tabulate the transcripts with alterantive first exons from the All_FilteredParsed output    
    '''
        
    AF = []    
    # Iterate through each transcript
    for transcript in df['transcript_id'].unique():
        class_transcript_exon,class_gencode_exon = class_by_transcript(transcript,  All_FilteredParsed)
        for i in class_gencode_exon:
            if "AF" in i:
                AF.append(transcript)
    
    return(AF)

