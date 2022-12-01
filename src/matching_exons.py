#!/usr/bin/env python
# Matching Exons

from .prepare_and_parse import class_by_transcript

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
    OtherClassLabels = ["ExtendedA5","ExtendedA3", "IR","IRMatch","TruncatedA5","TruncatedBothA3A5", "TruncatedA3"]
    
    # Iterate through each transcript to isolate the exon classifications
    AllKnownMatch = []
    SomeMatch = []
    Mis2Match = []
    OtherClass = []
    OtherClassMis2Match = []
    MisMatch = []
    MisjumpMatch = []
    MisjumpMatch_NotAll = []
    #for transcript in ["TALONT000983085"]:
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
        elif all(elem in gencode_exon for elem in Labels):
            if any(elem in OtherClassLabels for elem in class_split):
                if "MisMatch" in class_split or "MisjumpMatch" in class_split:
                    OtherClassMis2Match.append(transcript)
                else:              
                    OtherClass.append(transcript)
            elif "MisMatch" and "MisjumpMatch" in class_split:
                Mis2Match.append(transcript)
            elif "MisMatch" in class_split:
                MisMatch.append(transcript)
            elif "MisjumpMatch" in class_split and "MisMatch" not in class_split:
                MisjumpMatch.append(transcript)
            else:
                pass   
        elif not any(elem in OtherClassLabels for elem in class_split):
            MisjumpMatch_NotAll.append(transcript)        
        else:
            pass
    
    # Known transcripts from ONT reads identified using TALON
    TALON_identified_known_transcripts = list(filter(lambda x: ("ENS" in x), df["transcript_id"])) 
    AllKnownMatch = list(set(AllKnownMatch + TALON_identified_known_transcripts))
    
    print("Number of Transcripts with all exons and exact match:", len(AllKnownMatch))
    print("Number of Transcripts with all exons, with mismatch, misjump, exact match:", len(Mis2Match))
    print("Number of Transcripts with all exons, with mismatch and exact match:", len(MisMatch))
    print("Number of Transcripts with all exons, with misjumpmatch and exact match:", len(MisjumpMatch))
    print("Number of Transcripts with all exons, with other classification:", len(OtherClass))
    print("Number of Transcripts with all exons, but with mismatch and other classifications:", len(OtherClassMis2Match))
    print("Number of Transcripts with not all exons, but with match, mismatch, misjumpmatch:", len(MisjumpMatch_NotAll))
    print("Number of Transcripts with exact match but not all exons:", len(SomeMatch))
  
    return AllKnownMatch, Mis2Match, MisMatch, MisjumpMatch, OtherClass, OtherClassMis2Match, SomeMatch, MisjumpMatch_NotAll


def identify_all_exons(gencode, df, All_FilteredParsed):  
      
    '''
    Aim: Identify transcripts containing all exons independent of classifications 
    i.e. may have truncated A5, A3 or extended A5, A3 or intron retention
    :gencode = read gencode gtf
    '''
    
    # Generate the Labels for identification based on the gencode exon numbers
    max_exon = max([int(i) for i in gencode["updated_exon_number"]])
    Labels = ["Gencode_" + str(i) for i in range(1,max_exon+1)] 
    
    # Iterate through each transcript to isolate the exon classifications
    identified_transcripts = []
    for transcript in df['transcript_id'].unique():
        class_transcript_exon,class_gencode_exon = class_by_transcript(transcript, All_FilteredParsed)
        
        # remove "No" classifications
        class_gencode_exon = [s for s in class_gencode_exon if "No" not in s]
        gencode_num = [i.split("_",1)[1] for i in class_gencode_exon]
        
        # if equal or more exons than known exons and completely match all the labels (i.e all known exons)
        # then record
        if len(gencode_num) >= max_exon: 
            result =  all(elem in gencode_num for elem in Labels)
            if result:
                identified_transcripts.append(transcript)
    
    return(identified_transcripts)    
