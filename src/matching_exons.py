#!/usr/bin/env python
# Matching Exons

from prepare_and_parse import class_by_transcript

def identify_all_matching_exons(gencode, df, All_FilteredParsed):
    
    '''
    Aim: 
    1/ Identify transcripts with complete matching coordinates to known exons, must contain all exons
    2/ Identify transcripts with mismatching coordinates to known exons, must contain all exons
    :gencode = read gencode gtf
    '''
    
    # Generate the Matching Labels for identification based on the gencode exon numbers
    max_exon = max([int(i) for i in gencode["updated_exon_number"]])
    MatchLabels = ["Match_Gencode_" + str(i) for i in range(1,max_exon+1)]
    Labels = ["Gencode_" + str(i) for i in range(1,max_exon+1)] 
    OtherClassLabels = ["ExtendedA5","ExtendedA3","IR","IRMatch","TruncatedA5","TruncatedBothA3A5","TruncatedA3"]
    
    # Iterate through each transcript to isolate the exon classifications
    AllKnownMatch = []
    SomeMatch = []

    for transcript in df['transcript_id'].unique():
        class_transcript_exon,class_gencode_exon = class_by_transcript(transcript, All_FilteredParsed)
            
        # remove "No" classifications
        class_gencode_exon = [s for s in class_gencode_exon if "No" not in s]
        
        # set for parsing
        class_split = list(set([i.split("_",1)[0] for i in class_gencode_exon]))
        gencode_exon = list(set([i.split("_",1)[1] for i in class_gencode_exon]))
                  
        # if all the Match labels can be found in the classifications then record transcript id
        if all(elem in class_gencode_exon for elem in MatchLabels):
            AllKnownMatch.append(transcript)
        # All match to gencode exons, but do not contain all exons
        elif "Match" in class_split and all(ele == class_split[0] for ele in class_split):
            SomeMatch.append(transcript)
        else:
            pass
        
    print("Number of Transcripts with all exons and exact match:", len(AllKnownMatch))
    print("Number of Transcripts with exact match but not all exons:", len(SomeMatch))
  
    return AllKnownMatch, SomeMatch