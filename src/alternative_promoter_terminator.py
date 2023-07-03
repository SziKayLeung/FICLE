#!/usr/bin/env python
# alternative promoter and alternative terminator

import pandas as pd
from prepare_and_parse import class_by_transcript_pd

def identify_alternative_promoter(df, gencode, All_FilteredParsed):
    '''
    Transcripts with first exons that match first exons of shorter reference downstream are not also not considered AP 
    '''
    AP = []

    # Tabulate the corresponding exon in reference flattened structure of known first exons across all transcripts
    alt_1_exon = set(gencode.loc[gencode["altPromoters"],"updated_exon_number"].values)

    for transcript in df["transcript_id"].unique():
        # Record the corresponding gencode exon for the Transcript Exon 1
        #print(transcript)
        dat = class_by_transcript_pd(transcript, All_FilteredParsed)
        trans_exon1 = dat.loc[(dat["TranscriptExon"] == "Exon1") & (dat["Class"] != "No"),]
        # should only be 1 classification, unless there is in IR matched across
        if "IR" not in trans_exon1["Class"].values and "IRMatch" not in trans_exon1["Class"].values: assert len(trans_exon1) <= 1

        if len(trans_exon1) == 1:
            if str(trans_exon1["GencodeExon"].values[0]) in alt_1_exon: # if common values, then not AP
                AP.append(transcript)
        elif "IR" in trans_exon1["Class"].values:
            if str(min(trans_exon1["GencodeExon"].values)) in alt_1_exon:
                AP.append(transcript)
        else: # exon 1 not present
            AP.append(transcript)
    
    return AP
    
    
def identify_alternative_termination(df, gencode, All_FilteredParsed):
    
    alt_last_exon = set(gencode.loc[gencode["altTerminators"],"updated_exon_number"].values)

    AT = []
    for transcript in df["transcript_id"].unique(): 
        dat = class_by_transcript_pd(transcript, All_FilteredParsed)
        gencode_maxexon = max([int(i) for i in dat["GencodeExon"].values])

        if str(gencode_maxexon) not in alt_last_exon:
            AT.append(transcript)

    return AT
