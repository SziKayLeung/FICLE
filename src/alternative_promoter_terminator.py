#!/usr/bin/env python
# alternative promoter and alternative terminator

import pandas as pd
from prepare_and_parse import manual_gencode
from prepare_and_parse import class_by_transcript_pd

def identify_alternative_promoter(df, ES,gene, gencode, gencode_gtf, All_FilteredParsed):
    '''
    Aim: Identify alternative promoter based on the Exon Skipping table generated as the first exon identified
    :ES = Exon skipped table output from identify_exon_skipping() which also records position of first exon
    Note: Not alternative first exon which is subtle difference of the first exon
    '''
    output = []
    alt_clu = []
    for index,row in ES.iterrows():
        # if the first exon is not under Gencode Exon 1, then transcript characterised with alternative first exon
        if row["Gencode_1"] != "FirstExon" and row["Gencode_1"] != "IR":
            output.append(index)
    
    '''
    Transcripts with first exons that match first exons of shorter reference downstream are not also not considered AP 
    Exceptions
    Cd33 - Reference transcripts have first exon at "gencode 3" thus corresponding transcripts are not considered AP
    '''
    # list of coordinates and cooresponding exons of transcripts (not flattened)
    manualgencode_output = manual_gencode(gencode_gtf)
    Not_AP = []
    for transcript in df["transcript_id"].unique():
        # Record the corresonding gencode exon for the Transcript Exon 1
        dat = class_by_transcript_pd(transcript, All_FilteredParsed)
        trans_exon1 = [i.replace("Gencode_", "") for i in dat.loc[(dat["TranscriptExon"] == "Exon1") & (dat["Class"] != "No") ,"GencodeExon"].values]
        #print(trans_exon1)
        # Tabulate the corresponding exon in reference flattened structure of known first exons across all transcripts
        alt_1 = manualgencode_output[manualgencode_output["original_exon_number"] == 1]
        alt_1_exon = list(gencode.loc[(gencode["start"].isin(alt_1["start"].values)) & 
                         (gencode["end"].isin(alt_1["end"].values)),"updated_exon_number"].unique())
        #print(alt_1_exon)

        if set(trans_exon1).intersection(alt_1_exon): # if common values, then not AP
            Not_AP.append(transcript)

    if gene == "Cd33":
        if row["Gencode_3"] == "FirstExon": 
            alt_clu.append(index)  
    if gene == "Cd33": output = [x for x in output if x not in alt_clu]
    
    output = [x for x in output if x not in Not_AP]
    
    return(output)
    
    
def identify_alternative_termination(df, gencode, All_FilteredParsed):
    
    alt_lastexon = []
    for t in list(set(list(gencode["transcript_id"].values))):
        alt_lastexon.append(max(list(pd.to_numeric(gencode.loc[gencode["transcript_id"] == t,"updated_exon_number"]))))
        
    alt_lastexon = [str(i) for i in alt_lastexon]

    output = []
    #for transcript in ["PB.8560.365_TALONT001405949"]:
    for transcript in df["transcript_id"].unique(): 
        dat = class_by_transcript_pd(transcript, All_FilteredParsed)
        gencode_maxexon = max([int(i.replace("Gencode_"," ")) for i in dat["GencodeExon"].values])
        
        #print(gencode_maxexon)
        #print(alt_lastexon)

        if str(gencode_maxexon) not in alt_lastexon:
            output.append(transcript)

    return(output)
