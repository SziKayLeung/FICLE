#!/usr/bin/env python
# prepare_reference_transcript_gtf

from gtfparse import read_gtf
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import sys
import os

srcPath = os.path.dirname(os.path.realpath(__file__)) 
sys.path.insert(1, srcPath)

  
def determine_order(gencode):
    # extract only the first transcript of interest for determining order
    dat = gencode.loc[gencode["transcript_id"] == gencode["transcript_id"][1]]
    
     # sort by start coordinates in order depending of if sense or antisense 
    # note the exon_number are not in asending order from different transcripts
    if dat["start"][1] < dat["start"][len(dat)-1]: # sense i.e. first coordinate is smaller than last coordinate
        gencode = gencode.sort_values("start", ascending=True)
        order = "sense"
    elif dat["start"][1] > dat["start"][len(dat)-1]: # antisense
        gencode = gencode.sort_values("start", ascending=False)
        order = "antisense"
    else:
        print("Gencode incorrect order")
        
    return gencode, order 


def parse_gencode_reference(gencode_gtf, gene):
    
    '''
    Aim: Parse Gencode reference gtf of target gene and filter on exon genes
    *Important* Check this for every target gene that the final coordinates and exon number make sense
    
    1/ Obtain all coordinates of known exons 
    2/ Flatten gene structure by relabelling exon number based on the start and end coordinates of exons
    
    Scenario 1, 2, 3) Denoted the same exon with same start and/or end coordinates. 
    Scenario 4 also denoted the same exon, if within the start and end coordinates are within the reference exon 
    ref: ----exon------
    1:   ----exon------
    2:   ----exon--
    3:     --exon------
    4:     --exon--
    The reference exon refers to the exon before    
    
    '''
    
    # obtain coordinates of known exons, not exons without exon_number in genocode are dropped
    # only using transcripts that are protein-coding (removing lncRNA, antisense to gene etc)
    print("**** Extracting gtf", file=sys.stdout)
    print(gencode_gtf, file=sys.stdout)
    gencode = read_gtf(gencode_gtf)
    gencode = gencode[gencode["feature"] == "exon"]
    gencode = gencode[gencode["gene_type"] == "protein_coding"]
    gencode = gencode.loc[gencode["gene_name"] == gene,][["start","end","exon_number","transcript_id"]] 
    gencode.replace("", float("NaN"), inplace=True)
    gencode.dropna(subset = ["exon_number"], inplace=True)

    # sort the exons into order from start coordinates 
    # note: to avoid issues downstream whereby one exon could be annotated to different exon number of known transcripts 
    # essentially to flatten out all known transcripts into one

    # reindex after removing empty exon_number column
    gencode.index = np.arange(1, len(gencode) + 1)

    # sort by start coordinates in order depending of if sense or antisense 
    gencode, order = determine_order(gencode)

    # relabel the exon_number so the exon_number goes in assending order if the current value is < previous value in the column
    prev_final = gencode.loc[1,'exon_number'] 
    prev_start = int(gencode.loc[1,'start']) 
    prev_end = int(gencode.loc[1,'end']) 
    list_final = []
    for index, row in gencode.iterrows():
        exon = row["exon_number"]
        start = int(row["start"])
        end = int(row["end"])

        # if the current value is the same as the previous value
        # if the start or the end coordinate is the same as before, keep the previous value 
        # otherwise update
        if end == prev_end:
            list_final.append(prev_final)
        elif start == prev_start:
            list_final.append(prev_final)
        elif start > prev_start and end < prev_end and order == "sense":
            list_final.append(prev_final)
        elif start < prev_start and end > prev_end and order == "antisense":
            list_final.append(prev_final)
        else:
            prev_final = int(int(prev_final) + 1) 
            list_final.append(prev_final)

        prev_start = start
        prev_end = end

    # given creating an empty list, the first entry can be exon 1 or exon 2 
    # exon 1 if the top two entries belong to exon 1, or exon 2 if the second entry belongs to exon 2
    # therefore adjust the column in the later case by minus each element by 1, so starting from exon 1
    if int(list_final[0]) != 1:
        list_final[:] = [int(num) - 1 for num in list_final]

    gencode["updated_exon_number"] = list_final
    return(gencode)


def manual_gencode(gencode_gtf):
    gencode = read_gtf(gencode_gtf)
    gencode = gencode[gencode["feature"] == "exon"]
    gencode.index = np.arange(1,len(gencode)+1)
    
    output = []
    for transcript in gencode["transcript_id"].unique():
        df = gencode.loc[gencode["transcript_id"] == transcript,]
        count = 1
        for index,row in df.iterrows():
            output.append([row["start"],row["end"], count, transcript])
            count += 1
        #print(gencode.loc[gencode["transcript_id"] == transcript,][["start","end"]])
    
    return(pd.DataFrame(output, columns=["start","end","original_exon_number","transcript"]))


def replace_geneid(df, noISM):
    if "PB" in df["gene_id"][0]:
        print("Convert gtf file from SQANTI: PB gene id to gene name")
        
        # Extract the PB. gene id and create dictionary with the gene name from the classification file
        # Note this is from the subset classification file with 3'ISM removed
        # therefore will not replace PB.Gene ID that were removed, but these will not be useful for targeted data
        noISM["gene_id"] = ["PB." + i[1] for i in noISM["isoform"].str.split(".")]
        
        TargetGene = ["ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF"]
        dat = noISM[noISM["associated_gene"].isin(TargetGene)][["gene_id","associated_gene"]]
        dat.drop_duplicates(keep='first')
        gene_id_name = dict(zip(dat.gene_id,dat.associated_gene))
        
        # Replace PB gene id in the gtf file with the corresponding dictionary 
        df=df.replace({"gene_id": gene_id_name})
        return(df)
    else:
        return(df)
    

def parse_transcriptome_gtf(input_gtf, gene, noISM):
    
    '''
    Aim: Parse Input Gtf of target gene 
    '''
        
    # read gtf and filter only on exon (should only essentially remove "transcript" in gtf)
    df = read_gtf(input_gtf)
    
    # check if needs replacing the gene id column in the gtf from "PB.XX" to the gene name 
    df = replace_geneid(df, noISM)
    
    print("Original number of isoforms:", len(df["transcript_id"].unique()))
    df = df.loc[df["transcript_id"].isin(noISM["isoform"])]
    print("Retained number of isoforms after 3'ISM removal:", len(df["transcript_id"].unique()))
    
    # correct error in TALON gtf in misassigning transcripts for Apoe as ENSMUSG00000002981.10
    df.loc[(df["seqname"] == "chr7") & (df["start"] >= 19696109) & (df["end"] <= 19699188),"gene_id"] = "Apoe"
    
    # correct error in TALON gtf for having both Pt2kB and Pt2kb classifications
    df.loc[df["gene_id"] == "Ptk2B","gene_id"] = "Ptk2b"
    
    print("**** Extracting for transcripts associated with:", gene)
    df = df[df["gene_id"] == gene]
    df = df[df["feature"] == "exon"]
    print("Number of Detected transcripts:", len(df["transcript_id"].unique()))
   
    return(df)

   
def parse_transcript(gencode, parsed_transcriptome, transcript, wobble, IR_threshold):

    '''
    Aim: parse through each transcript to classify the splicing events based on coordinates
    : gencode = parsed gencode using the updated_exon_number column 
    : transcript = the transcript ID 
    : wobble = the maximum number of bp from the splice junction (SJ) to still be called as the same 
    : IR_threshold = the maximum number of bp from the SJ to be classified as intron retained 
    
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
    # Target end coordinates exact but Target start coordinates less than wobble range but more than IR range ==> "TruncatedA5"
    # Target start and exon coordinates less than the wobble range ==> "TruncatedBothA3A5"
    
    ** Extended **
    # Target start coordintes exact, but end coordinates more than wobble range but less IR range ==> "ExtendedA3"
    # Target end coordinates exact, but start coordinates more than wobble range but less than IR range ==> "ExtendedA5"
    
    ** Intron Retention **
    # Target start coordinates exact, but end coordinates more than IR range allowed ==> IR
    # Target end coordinates exact, but start coordinates more than IR range ==> IR
    # Start end coordinates exact, but end coordinates match other downstream gencode end ==> IRMatch
    
    ** None **
    # No match of the start or end coordinates to the reference gencode exon
    '''
    
    # info to extract from the gencode reference 
    # maximum exon number and all the ends of exons
    max_exon = max([int(i) for i in gencode["updated_exon_number"]])
    gencode["updated_exon_number"] = [str(i) for i in gencode["updated_exon_number"]]
    #print("Number of exons known:", max_exon)
    list_gencode_ends = gencode[["end"]]
    
    # Filter 
    parsed = []  
    filtered = parsed_transcriptome[parsed_transcriptome["transcript_id"] == transcript]
    
    # determine if gene is sense or antisense
    gencode, order = determine_order(gencode)
    
    if order == "sense": # sense i.e. first coordinate is smaller than last coordinate
        direction = "forward"
    else:
        direction = "reverse"
        filtered = filtered.iloc[::-1] # reverse the rows
    
    # Reindex
    filtered.index = np.arange(1, len(filtered) + 1)
    
    # Iterature through each exon
    for index, row in filtered.iterrows():
        
        # capturing transcript start and end coordinates
        transcript = row["transcript_id"]
        transcript_START = int(row["start"])
        transcript_END = int(row["end"])
        transcript_EXON = str(index)
        
        # Match each exon coordinate of detected target transcripts to each exon in reference gencode
        for gencode_index, gencode_row in gencode.iterrows():
            
            # capturing gencode reference start and end coordinates 
            gencode_START = int(gencode_row["start"])
            gencode_END = int(gencode_row["end"])
            gencode_EXON = str(gencode["updated_exon_number"][gencode_index])
            
            # For Intron retention classification 
            # Find the gencode end sites that have matching start sites to the target transcript start
            paired_ends = list(gencode[gencode["start"] == transcript_START]["end"])
            # Identify the remaining gencode end sites for further intron retention classification
            other_gencode_ends = list(set(list_gencode_ends["end"]) - set(paired_ends))
            
            # Exact Coordinates ==> MATCH    
            if(transcript_START == gencode_START and transcript_END == gencode_END):
                parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";Match_Gencode_" + gencode_EXON)
                
            # Exon Coordinates within wobble range ==> MATCH
            elif (gencode_END - wobble <= transcript_END <= gencode_END + wobble and 
                  gencode_START - wobble <= transcript_START <= gencode_START + wobble):
                parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";Match_Gencode_" + gencode_EXON)
            
            # Start coordinates the same, but transcript end shorter than gencode end ==> TruncatdA3
            # but not as far as the gencode end with the IR threshold
            #elif (transcript_START == gencode_START and gencode_END - IR_threshold <= transcript_END <= gencode_END - wobble):
            elif (gencode_START - wobble <= transcript_START <= gencode_START + wobble and transcript_END <= gencode_END - wobble):
                if direction == "forward": parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";TruncatedA3_Gencode_" + gencode_EXON)
                if direction == "reverse": parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";TruncatedA5_Gencode_" + gencode_EXON)
            
            # End coordinates the same, but transcript start shorter than gencode start ==> TruncatedA5
            # but not as far as the gencode start with the IR threshold
            #elif (transcript_END == gencode_END and gencode_START + IR_threshold >= transcript_START >= gencode_START - wobble):
            elif (gencode_END - wobble <= transcript_END <= gencode_END + wobble and transcript_START >= gencode_START - wobble):
                if direction == "forward": parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";TruncatedA5_Gencode_" + gencode_EXON) 
                if direction == "reverse": parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";TruncatedA3_Gencode_" + gencode_EXON)
            
            # End and Start Exon coordinates shrunk within the gencode start and end exon coordinates ==> TruncatedBothA3A5
            elif (transcript_END <= gencode_END - wobble and gencode_END >= transcript_START > gencode_START - wobble):
                parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";TruncatedBothA3A5_Gencode_" + gencode_EXON)
            
            elif (gencode_END + wobble <= transcript_END <= gencode_END + IR_threshold and 
                  gencode_START - IR_threshold <= transcript_START <= gencode_START - wobble or
                  gencode_END - wobble >= transcript_END >= gencode_END - IR_threshold and
                  gencode_START + IR_threshold >= transcript_START >= gencode_START + wobble):
                parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";ExtendedBothA3A5_Gencode_" + gencode_EXON)
            #elif (transcript_END == gencode_END and gencode_row["exon_number"] == 1):
            #    parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";UTR_Gencode_" + gencode_EXON)
            
            # Start coordinates the same, but end coordinates more than the wobble range allows but less than the IR 
            # thus not considered intron retention
            elif (gencode_START - wobble <= transcript_START <= gencode_START + wobble and gencode_END + wobble <= transcript_END <= gencode_END + IR_threshold):
                if direction == "forward": parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";ExtendedA3_Gencode_" + gencode_EXON)
                if direction == "reverse": parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";ExtendedA5_Gencode_" + gencode_EXON)
            
            # End coordinates the same, but start coordinates more than the wobble range allows but less than the IR 
            # thus not considered intron retention
            elif (gencode_START - IR_threshold <= transcript_START <= gencode_START - wobble and gencode_END - wobble <= transcript_END <= gencode_END + wobble): 
                if direction == "forward": parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";ExtendedA5_Gencode_" + gencode_EXON)
                if direction == "reverse": parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";ExtendedA3_Gencode_" + gencode_EXON)
            
            # Start site the same but the end coordinates match downstream exon coordinates rather than the corresponding gencode end
            elif (transcript_START == gencode_START and transcript_END in other_gencode_ends):
                IR_exon = gencode[gencode["end"] == transcript_END]
                parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON + ";IRMatch_Gencode_" + gencode_EXON)
                parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON + ";IRMatch_Gencode_" + max(IR_exon["updated_exon_number"]))
            
            # Start coordinates the same with wobble but end coordinates further than the IR threshold allowed OR
            # End coordinates the same with wobble but start coordinates further starts than the IR threshold allowed                
            elif (gencode_START - wobble <= transcript_START <= gencode_START + wobble and 
                  transcript_END >= gencode_END + IR_threshold or
                  gencode_START - IR_threshold >= transcript_START and 
                  gencode_END - wobble <= transcript_END <= gencode_END + wobble):
                parsed.append(transcript + ";Transcript_Exon"+ transcript_EXON +";IR_Gencode_" + gencode_EXON)

            else: # No classification to the reference gencode exon
                parsed.append(transcript +  ";Transcript_Exon"+ transcript_EXON +";No_Gencode_" + gencode["updated_exon_number"][gencode_index])
              
    return(parsed)



def filter_parsed_transcript(gene, gencode, parsed, output_log):

    '''
    Aim: Filter the parsed output to prioritise the output 
    Note: due to multiple reference transcripts with same or similar exons, 
    the parse function generates redundant classifications
    i.e one target exon can be complete match to one reference exon but truncated/extended to another reference exon
    :gencode: parsed output gencode gtf from parse_gencode_reference()
    :parsed: output from parse_transcript()

    
    On the parsed output i.e. for each transcript,
    1/ Iterate through each gencode exon matched and identify the classifications that correspond to it
    2/ If more than one classification for each gencode exon, then filter 
    
    Parsed element structure: 
    <transcript_id>;<transcript_exon>;<class_gencode>
    <transcript_id>;Transcript_Exon<ExonNum>;<exon_class>_Gencode>_<ExonNum>
    
    *** Filter Priorities ***
    # If there is a match for that gencode exon, remove all other classifications 
    # If there is a "No match" for that gencode exon, but there are transcript exons that do match 
        # i.e other classifications, remove the "No" elements
    # If there is TruncatedBothA3A5 for that exon in the presence of TruncatedA3 or TruncatedA5 
        # then remove TruncatedBothA3A5 as the exact match from TruncatedA3 and TruncatedA5 more important  
    # If there is IR for that exon in the presence of TruncatedA3/TruncatedA5/ExtendedA3/A5
        # them remove IR, as the distance from truncation/extension is smaller than the IR threshold
        # suggesting that there are other exons from other transcripts that are more similar to exon of interest
    # If there is an IR Match, then the start coordinates most important and so overrides ExtendedA5
    '''
    
    # info to extract from the gencode reference 
    # maximum exon number to generate the list of possible results from the parsed output
    # info to extract from the gencode reference 
    # maximum exon number to generate the list of possible results from the parsed output
    max_exon = max([int(i) for i in gencode["updated_exon_number"]])
    Gencode_exons = ["Gencode_" + str(i) for i in range(1,int(max_exon) + 1)]
    set_class = pd.DataFrame([[i.split(";",3)[0],
                                   i.split(";",3)[1].split("_",2)[1], 
                                   i.split(";",3)[2].split("_",2)[0], 
                                   i.split(";",3)[2].split("_",1)[1]] for i in parsed], 
                             columns=["Transcript","TranscriptExon", "Class","GencodeExon"])


    remove_class = [] # list of parse elements to filter out

    # iterate through each gencode exon to remove redundant classifications to that exon 
    for g in Gencode_exons:
            #print(g)
            redundant = set_class[(set_class['GencodeExon']==g)]
            
            # if multiple classifications for the gencode exon
            if len(redundant) >= 1:
                # work on that specific Gencode exon and further split by class
                exon_class_redundant = list(redundant["Class"])
                #print(exon_class_redundant)

                # if there is a "Match" ==> remove all other classifications to that gencode exon
                if "Match" in exon_class_redundant:
                    remove = list(filter(lambda x: not "Match" in x, exon_class_redundant)) 
                    remove_class.extend([i + "_" + g for i in list(filter(lambda x: not "Match" in x, exon_class_redundant))])

                # if there a "No" but there are other classifications for that same gencode exon
                # ==> remove the no classifications 
                if "No" in exon_class_redundant:
                    if exon_class_redundant.count("No") != len(exon_class_redundant):
                        remove_class.extend([i + "_" + g for i in list(filter(lambda x: "No" in x, exon_class_redundant))])
                
                if "IRMatch" in exon_class_redundant:
                    remove_class.extend([i + "_" + g for i in list(filter(lambda x: "ExtendedA3" in x, exon_class_redundant))])
                    remove_class.extend([i + "_" + g for i in list(filter(lambda x: "ExtendedA5" in x, exon_class_redundant))])
                    remove_class.extend([i + "_" + g for i in list(filter(lambda x: "IR" in x, exon_class_redundant))])
                    remove_class = [i for i in remove_class if "IRMatch" not in i]

                # if there is a Truncated both A3A5 and TruncatedA3 or TruncatedA5 then remove
                if "TruncatedA3" in exon_class_redundant or "TruncatedA5" in exon_class_redundant:
                    if "TruncatedBothA3A5" in exon_class_redundant:
                #        remove_class.extend([i + "_" + g for i in list(filter(lambda x: "TruncatedBothA3A5" in x, exon_class_redundant))])
                        remove_class.extend([i + "_" + g for i in list(filter(lambda x: "IR" in x, exon_class_redundant))])
                

                if "ExtendedA3" in exon_class_redundant or "ExtendedA5" in exon_class_redundant:
                        remove_class.extend([i + "_" + g for i in list(filter(lambda x: "IR" in x, exon_class_redundant))])                           

        
                # Further filtering 
                g_df = set_class.loc[set_class["GencodeExon"] == g,]        
                #print(set(g_df["Class"]))
                remove_class_lst = ["TruncatedBothA3A5_" + g, "TruncatedA5_" + g,"TruncatedA3_" + g, "ExtendedA3_" + g, "ExtendedA5_" + g, "IR_" + g]
                remove_class_lst_2 = ["TruncatedBothA3A5_" + g, "TruncatedA5_" + g,"TruncatedA3_" + g, "ExtendedA3_" + g, "ExtendedA5_" + g]


                if "TruncatedA5" in list(g_df["Class"]) and "TruncatedA3" in list(g_df["Class"]):
                    A3 = [int(i.replace("Exon","")) for i in list(g_df.loc[g_df["Class"] == "TruncatedA3","TranscriptExon"])]
                    A5 = [int(i.replace("Exon","")) for i in list(g_df.loc[g_df["Class"] == "TruncatedA5","TranscriptExon"])]

                    if set(A3).intersection(A5):
                        com = list(set(A3).intersection(A5))
                        remove_class.extend(remove_class_lst) 
                        parsed.append(set_class["Transcript"][0] + ";Transcript_Exon" + str(com[0]) + ";MisMatch_" + g)

                    #if len(set(A3)) and len(set(A5)) == 1:
                    #    if abs(A3[0] - A5[0]) == 1:
                    #        remove_class.extend(remove_class_lst) 
                    #        parsed.append(set_class["Transcript"][0] + ";Transcript_Exon" + str(A3[0]) + ";MisjumpMatch_" + g)
                    #        parsed.append(set_class["Transcript"][0] + ";Transcript_Exon" + str(A5[0]) + ";MisjumpMatch_" + g)
                
                
                if "TruncatedA5" in list(g_df["Class"]) and "TruncatedBothA3A5" in list(g_df["Class"]):
                    A5 = [int(i.replace("Exon","")) for i in list(g_df.loc[g_df["Class"] == "TruncatedA5","TranscriptExon"])]
                    A3A5 = [int(i.replace("Exon","")) for i in list(g_df.loc[g_df["Class"] == "TruncatedBothA3A5","TranscriptExon"])]
                    #if len(set(A5)) and len(set(A3A5)) == 1:
                    #    if abs(A5[0] - A3A5[0]) == 1:
                    #        remove_class.extend(remove_class_lst) 
                    #        parsed.append(set_class["Transcript"][0] + ";Transcript_Exon" + str(A5[0]) + ";MisjumpMatch_" + g)
                    #        parsed.append(set_class["Transcript"][0] + ";Transcript_Exon" + str(A3A5[0]) + ";MisjumpMatch_" + g)
                
                if g == "Gencode_1":
                    if "Match" or "MisjumpMatch" or "MisMatch" or "IRMatch" in exon_class_redundant:
                        pass
                    elif exon_class_redundant.count("No") == len(exon_class_redundant):
                        pass
                    else:
                        # remove all classes except for IR for first Exon as still need to tabulate 
                        remove_class.extend(remove_class_lst_2)
                        min_transcript_exon = min(list(i.replace("Exon","") for i in g_df["TranscriptExon"]))
                        parsed.append(set_class["Transcript"][0] + ";Transcript_Exon" + min_transcript_exon + ";AF_" + g)


    '''
    Aim: Final step of filtering - each transcript exon only corresponds to one gencode exon 
    1st priority is placed on exons that have the same number i.e Transcript exon 1 and Gencode exon 1 
    1/ Create a dataframe from the filtered list 
    2/ Loop through each transcript exon, and > 1 entry (i.e. duplicated records), then select the row with the same gencode exon
    **Note** if the transcript exon number is not within the gencode exon number, then pass for exception of transcripts with novel exons beyond first
    ==> Different exon order to the reference even if matching downstream
    3/ Remove other duplicated records


    # Filter most important for removing transcript exons in these situations:
    Transcript:                                [Exon]
    Gencode reference A: [Exon]--------------[  Exon]
    Gencode reference B: [Exon                      ]
    The transcript exon would correspond to before gencode reference but different exon number i.e. Exon 2 with ReferenceA than Exon 1 with ReferenceB
    
    **Note**
    Transcript:               [Exon]---[Exon]------[1st Exon]-----[Exon]
    Gencode reference:                             [1st Exon]-----[Exon]
    '''

    # convert list to dataframe 
    parsed_df = pd.DataFrame([[i, 
                                 i.split(";",3)[0],
                                 i.split(";",3)[1].split("_",2)[1], 
                                 i.split(";",3)[2].split("_",2)[0], 
                                 i.split(";",3)[2].split("_",1)[1],
                                 i.split(";",1)[1]] for i in parsed], 
                             columns=["Original","Transcript","TranscriptExon", "Class","GencodeExon","Labels"])

     # loop through each exon in the dataframe
    Ext_Tru_Class = ["TruncatedA3", "TruncatedA5", "ExtendedA3","ExtendedA5"]
    for exon in set(parsed_df["TranscriptExon"]):
        # if more than one record
        if len(parsed_df.loc[parsed_df["TranscriptExon"] == exon]) > 1:
            redundant_gencode_class = list(parsed_df.loc[parsed_df["TranscriptExon"] == exon,"Class"]) 
            
            #print(redundant_gencode_class)
            # If there is a match, then remove all other records
            if "Match" in redundant_gencode_class:
                remove_class.extend(list(parsed_df.loc[(parsed_df["TranscriptExon"] == exon) & (parsed_df["Class"] != "Match")]["Labels"]))
        
            # remove No classifications if other classifications are present
            if parsed_df['Class'].value_counts().No != len(parsed_df):
                remove_class.extend(list(parsed_df.loc[(parsed_df["TranscriptExon"] == exon) & (parsed_df["Class"] == "No")]["Labels"]))
            
            # Extension/Truncation with one matched end prioritised over truncated both ends 
            if any(elem in Ext_Tru_Class for elem in redundant_gencode_class):
                remove_class.extend(list(parsed_df.loc[(parsed_df["TranscriptExon"] == exon) & (parsed_df["Class"] == "TruncatedBothA3A5")]["Labels"]))
                remove_class.extend(list(parsed_df.loc[(parsed_df["TranscriptExon"] == exon) & (parsed_df["Class"] == "ExtendedBothA3A5")]["Labels"]))
                if gene != "Trem2":
                    remove_class.extend(list(parsed_df.loc[(parsed_df["TranscriptExon"] == exon) & (parsed_df["Class"] == "IR")]["Labels"]))   
            
            if "IRMatch" in redundant_gencode_class:
                remove_class.extend(list(parsed_df.loc[(parsed_df["TranscriptExon"] == exon) & (parsed_df["Class"] == "TruncatedBothA3A5")]["Labels"]))
                remove_class.extend(list(parsed_df.loc[(parsed_df["TranscriptExon"] == exon) & (parsed_df["Class"] == "TruncatedA3")]["Labels"]))
                remove_class.extend(list(parsed_df.loc[(parsed_df["TranscriptExon"] == exon) & (parsed_df["Class"] == "TruncatedA5")]["Labels"]))
                
            if "MisjumpMatch" in redundant_gencode_class:
                remove_class.extend(list(parsed_df.loc[(parsed_df["TranscriptExon"] == exon) & (parsed_df["Class"].isin([Ext_Tru_Class]))]["Labels"]))

            if "TruncatedBothA3A5" in redundant_gencode_class:
                remove_class.extend(list(parsed_df.loc[(parsed_df["TranscriptExon"] == exon) & (parsed_df["Class"] == "IR")]["Labels"])) 
    
    
    # remove labels within the dataframe and convert back to list    
    parsed_df = parsed_df[~parsed_df["Labels"].isin(remove_class)]
    parsed_df = parsed_df.drop_duplicates()
    
    # remove duplicate exon in transcript 
    for exon in set(parsed_df["TranscriptExon"]):
                exon_num = exon.replace("Exon","") 
            
                # record the gencode exon
                redundant_gencode_exon = parsed_df.loc[parsed_df["TranscriptExon"] == exon,"GencodeExon"].values
                redundant_gencode_exon = [i.replace("Gencode_","") for i in redundant_gencode_exon] 
                
                # only if the exon_number is within the list of gencode_exon number 
                if exon_num in redundant_gencode_exon: 
                    if len(redundant_gencode_exon) > 1: 
                        for i in redundant_gencode_exon:
                            # the X from Exon_X and Gencode_X are the same, then keep that class, and remove other classes 
                            if exon_num != i:
                                specific_class = list(parsed_df[(parsed_df["TranscriptExon"] == exon) & (parsed_df["GencodeExon"] == "Gencode_" + i)]["Class"])
                                for c in specific_class:
                                    remove_class.extend(["Transcript_" + exon + ";" + c + "_Gencode_" + i])

    parsed_df = parsed_df[~parsed_df["Labels"].isin(remove_class)]
    
    # To account for complexity of Tardbp 
    if gene != "Tardbp":    
        for exon in set(parsed_df["GencodeExon"]):
                exon_num = exon.replace("Gencode_","") 

                # record the transcript exon
                redundant_transcript_exon = parsed_df.loc[parsed_df["GencodeExon"] == exon,"TranscriptExon"].values
                redundant_transcript_exon = [i.replace("Exon","") for i in redundant_transcript_exon] 
                #print(redundant_transcript_exon)

                # only if the exon_number is within the list of gencode_exon number 
                if exon_num in redundant_transcript_exon: 
                    if len(redundant_transcript_exon) > 1: 
                        for i in redundant_transcript_exon:
                            # the X from Exon_X and Gencode_X are the same, then keep that class, and remove other classes 
                            if exon_num != i: 
                                specific_class = list(parsed_df[(parsed_df["GencodeExon"] == exon) & (parsed_df["TranscriptExon"] == "Exon" + i)]["Class"])
                                for c in specific_class:
                                    remove_class.extend(["Transcript_Exon" + i + ";" + c + "_Gencode_" + exon_num])


        parsed_df = parsed_df[~parsed_df["Labels"].isin(remove_class)]
        parsed = parsed_df['Original'].to_list()  
    
    # Iterate through the parsed list and remove other remove_class elements
    filtered = []
    for s in parsed:
        s_split = s.split(';',3)[2] 
        if s_split not in remove_class:
            filtered.append(s)
    
    # write the output to a log file 
    for element in filtered:
        output_log.write(element + "\n")  

    filtered = list(set(filtered))
    return(filtered)
#filter_parsed_transcript(gene, gencode, parse_transcript(gencode, df, "TALONT001069788", 10, 100), output_log)


def class_by_transcript(transcript, All_FilteredParsed):
    
    '''
    Aim: Quickly filter output from filter_parsed_transcript() using the transcript ID
    # Used in multiple functions downstream
    :transcript = transcript ID
    
    Output: 
    :class_transcript_exon = list of transcript exons classified 
    :class_gencode_exon = list of corresponding gencode exons classified with the classification types
    '''
    
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
                                 i.split(";",3)[0],
                                 i.split(";",3)[1].split("_",2)[1], 
                                 i.split(";",3)[2].split("_",2)[0], 
                                 i.split(";",3)[2].split("_",0)[0],
                                 i.split(";",3)[2].split("_",1)[1],
                                 i.split(";",1)[1]] for i in All_FilteredParsed], 
                             columns=["Original","Transcript","TranscriptExon", "Class","ClassGencodeExon","GencodeExon","Labels"])
    
    output_df = parsed_df.loc[parsed_df["Transcript"] == transcript,]
    return(output_df)


def generate_split_table(list2table,colname):
    
    # Create as dataframe, split the column by comma into two 
    df = pd.DataFrame({'Cat':list2table})
    df[['transcript_id', colname]] = df['Cat'].str.split(',', 1, expand=True)
    
    return(df)


def reidentify_isoform_dataset(input_file):
    
    '''
    #### Function only applicable for merged datasets ###
    Aim: Define from the isoforms of a merged dataset (i.e. PacBio and ONT) from the respective dataset
    params:
        input_file = dataframe contains an "isoform" column with the PacBio and ONT isoforms merged 
        ** Important **: PacBio ID always before ONT ID if isoform present in both datasets, separated by "_"
    
    Output: Appends two additional columns (dataset, and ONT isoforms)
    '''

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
