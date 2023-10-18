#!/usr/bin/env python
# A5A3 and AF

import pandas as pd
from prepare_and_parse import class_by_transcript
from prepare_and_parse import class_by_transcript_pd
from prepare_and_parse import generate_split_table

def identify_A5A3_transcripts(args, df, All_FilteredParsed):
    
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
        output_df = generate_split_table(output_lst,"A5A3")
        output_transcripts = list(output_df["transcriptID"].unique())

        # Output Stats
        output_df[['category', 'gencodeExon']] = output_df['A5A3'].str.split('_', 1, expand=True)
        output_df = output_df[['transcriptID','category','gencodeExon']]
        output_df.to_csv(args.gene_stats_dir + args.genename + "_A5A3_transcript_level.csv",index = False)

        output_pertrans_counts = output_df.groupby(['transcriptID']).size().to_frame('numEvents')
        output_pertrans_counts.to_csv(args.gene_stats_dir + args.genename + "_A5A3_transcript_counts.csv")

        output_counts = output_df.groupby("category").size().to_frame('numTranscripts')
        output_counts.to_csv(args.gene_stats_dir + args.genename + "_A5A3_events_counts.csv")

    except:
        print("No Transcripts with A5A3 truncation or extension")
    return output_df, output_transcripts

    
def identify_alternative_first(df, All_FilteredParsed):
    
    '''
    Aim: Tabulate the transcripts with alterantive first exon (extended,truncated) from the All_FilteredParsed output
    
    Output: 
    output_df: Table of Transcripts with corresponding first exon with extension or truncation
    output_counts: The total number of transcripts with extended or truncated first exons
    output_transcripts: List of transcripts with AF
    '''
        
    AF = []    
    for transcript in df['transcript_id'].unique():
        dat = class_by_transcript_pd(transcript,  All_FilteredParsed)
        datAF = dat[(dat['Class'].str.contains('Truncated|Extended', case=False, regex=True)) & (dat['TranscriptExon'] == 'Exon1')]
        if len(datAF) > 0:
            AF.append([transcript + "," + datAF["ClassGencodeExon"].values[0]])
    
    output_df = pd.DataFrame()
    output_counts = pd.DataFrame()
    output_transcripts = []
    try:
        output_df = generate_split_table(sum(AF, []),"AF")[["transcriptID","AF"]]
        output_counts = output_df['AF'].str.split('_', 1, expand=True).groupby([0]).size().reset_index(name='count')
        output_transcripts = list(output_df["transcriptID"].unique())     
    except:
        print("No Transcripts with AF")
    
    return output_df, output_counts, output_transcripts

