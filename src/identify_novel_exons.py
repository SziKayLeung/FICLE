#!/usr/bin/env python
# alternative promoter and alternative terminator

import pandas as pd
import numpy as np
from prepare_and_parse import class_by_transcript_pd
from prepare_and_parse import generate_split_table
from prepare_and_parse import determine_order

def novel_exon_stats(args, NE, NE_classify):
    
    '''
    Aim: Tabulate Stats for Novel Exons
    :NE = Novel exon table output from identify_novel_exon()
    :NE_classify = Further classification of novel exon table output from classify_novel_exon()
    
    Output: 
    NE_pertrans_counts = Number of novel exons per transcript
    NE_classify_counts = The total number of transcripts with beyond first, internal and beyond last
    NE_pertrans_classify_counts = Number of novel exons (classified by internal etc) per transcript
    '''
    
    NE_pertrans_classify_counts = pd.DataFrame()
    if len(NE) > 0:
        # Table 1: Number of novel exons per transcript
        # From the NE table, group by the number of transcripts
        NE_pertrans_counts = NE.groupby(['transcriptID'])['novelexon'].count().reset_index()
        NE_pertrans_counts.columns = ["transcriptID","numNovelExons"]
        NE_pertrans_counts.to_csv(args.gene_stats_dir + args.genename + "_NE_transcript_counts.csv", index = False)

        # Output 3: 
        NE_classify_counts = NE_classify.groupby(['classification']).count().reset_index()
        NE_classify_counts.columns = ["typeNovelExon","numTranscripts"]
        NE_classify_counts.to_csv(args.gene_stats_dir + args.genename + "_NE_type_counts.csv", index = False)

        # Output 4:
        NE_pertrans_classify_counts = NE_classify.groupby(['classification','transcriptID']).size().reset_index()
        NE_pertrans_classify_counts.columns = ['NEtype', 'transcriptId', 'numNovelExons']
        NE_pertrans_counts.to_csv(args.gene_stats_dir + args.genename + "_NE_transcript_counts.csv", index = False)
        print("Total Number of transcripts with novel exon: ", len(NE_pertrans_counts["transcriptID"]))
    
    return NE_pertrans_classify_counts


def identify_novel_exon(args, df, gencode, All_FilteredParsed):
        
    
    '''
    Aim: 
    1/ Identify transcripts with novel exons by identifying the exon that has no classifications in the All_FilteredParsed
    2/ Identify coordinates of novel exons i.e. number of novel exons since multiple transcripts can have the same novel exon
    :df = read transcriptome gtf    
    '''

    # Record the total number of exons in each detected transcript for downstream 
    detected_trans_exon_counts = df.groupby(['transcript_id'])['feature'].count().reset_index()
    detected_trans_exon_counts_dict = dict(zip(list(detected_trans_exon_counts["transcript_id"]),list(detected_trans_exon_counts["feature"])))

    # All possible classifications
    AllLabels = ["AF", "Match", "MisMatch", "MisjumpMatch", "ExtendedA5","ExtendedA3",
                 "UTR_Gencode", "FirstExon","TruncatedA5","TruncatedBothA3A5", "TruncatedA3","IR","IRMatch"]
    max_exon = max([int(i) for i in gencode["updated_exon_number"]])
    Gencode_exons = ["Gencode_" + str(i) for i in range(1,max_exon+1)] 

    output = []
    #for transcript in ["TALONT000703466"]:
    for transcript in df['transcript_id'].unique():
        detected_max_exons = [val for key, val in detected_trans_exon_counts_dict.items() if transcript in key][0]
        trans_pd = class_by_transcript_pd(transcript, All_FilteredParsed)

        for exon in range(1,detected_max_exons + 1): # + 1 to make sure the final exon also captured 
            Gencode_split = list(trans_pd.loc[trans_pd["TranscriptExon"] == "Exon" + str(exon),"Class"].values)

            result =  any(elem in AllLabels  for elem in Gencode_split)
            if result:
                pass
            else:
                # include chr, start and end coordinates 
                # exon - 1 to account index from 0 in df
                seqname = str(df.loc[df["transcript_id"] == transcript, "seqname"].values[exon - 1])
                startCo = str(df.loc[df["transcript_id"] == transcript, "start"].values[exon - 1])
                endCo = str(df.loc[df["transcript_id"] == transcript, "end"].values[exon - 1])
                output.append(transcript + "," + str(exon) + "," + seqname + ":" + startCo + "-" + endCo + "," + startCo)

    output_df = pd.DataFrame()
    NE_coordinates = []

    try:
        output_df = pd.DataFrame([line.split(',') for line in output], columns=['transcriptID', 'novelexon', 'co-ordinates', 'start'])
        # sort by ascending start coordinates
        output_df = output_df.sort_values(by=['start']) 
        output_df = output_df[["transcriptID","novelexon","co-ordinates"]]
        print("Number of unique novel exons:" , len(set(output_df["co-ordinates"]))) 
    
    except:
            print("no novel exons")
    
    return output_df


def classify_novel_exon(args, gencode, order, df, NE, All_FilteredParsed):    
    
    '''
    Aim: Further classify location of novel exon as being beyond the first exon, internal or beyond the last exon
    # Beyond First: Novel exon is to the left of the first known exon i.e. START exon < Minimum known START exon 
    # Internal: Novel exon is between the known gene coordinates 
    # Beyond Last: Novel exon is to the right of the last known exon i.e. START exon > Maximum known END exon
    # First: Novel exon is modification of the first exon but with coordinates that fall out of A5'A3'
    # Last: Novel exon is the final exon of the transcript but not beyond last known genome exon
    :df = read transcriptome gtf
    :NE = novel exon table output from identify_novel_exon()
    
    1/ Find the coordinates of the novel exon associated to each transcript 
    2/ Compare the coordinates to the known coordinates for further classification
    '''
    # Gencode start and end exons for determining location of novel exon    
    if order == "antisense":
        df = df.sort_values("start", ascending=False)
        gene_start = max(gencode["start"])
        gene_end = min(gencode["end"])
    else:
        gene_start = min(gencode["start"])
        gene_end = max(gencode["end"])

    # empty dataframe and list
    output_df = pd.DataFrame()
    NExons_BeyondFirst = []
    NExons_BeyondFirstLast = []
    NExons_BeyondLast = []
    NExons_Internal = []
    NExons_First = []
    NExons_Last = []

    if len(NE) > 0:
        # Iterate through each transcript detected to have novel exon 
        # Subset the transcriptome with coordinates corresponding to that transcript
        # Re-index to find the coordinates of the novel exon and compare to the known gencode START and END

        # Beyond_First; Beyond_First_Last, Beyond_Last; Internal_NovelExon; First 
        NovelExons_Subcategory = []
        #NE = NE.loc[NE["transcript_id"] == "PB.22007.5903"]
        for index,row in NE.iterrows():
            transcript = row["transcriptID"]
            # subset the transcriptome and reset index to find coordinates of novel exon
            subset = df[df["transcript_id"] == transcript].reset_index()
            # strip the column from the novel exon table for the number 
            novel_exon = str(row["novelexon"])
            # -1 to consider the index starts at 0
            novel_exon_start = subset.loc[int(novel_exon) -1,"start"]
            novel_exon_end = subset.loc[int(novel_exon) -1 ,"end"]
            
            if order == "sense":
                if novel_exon_end < gene_start:
                    NovelExons_Subcategory.append(transcript + ",Beyond_First")
                elif novel_exon_start > gene_start and novel_exon_end < gene_end:
                    NovelExons_Subcategory.append(transcript + ",Internal_NovelExon")
                elif novel_exon_start > gene_end:
                    NovelExons_Subcategory.append(transcript + ",Beyond_Last")
                elif novel_exon == "1":
                    NovelExons_Subcategory.append(transcript + ",First")
                else:
                    # check if the novel exon is the last exon (
                    # i.e. missing exon is beyond the final known exon in the transcript)
                    trans_pd = class_by_transcript_pd(transcript, All_FilteredParsed)
                    maxexondetected = max([int(i.replace("Exon","")) for i in trans_pd["TranscriptExon"].values])
                    if int(novel_exon) > maxexondetected:
                        NovelExons_Subcategory.append(transcript + ",Last")
                    else:
                        print("Transcript not classified for Novel Exons:", transcript)
            else:
                if novel_exon_end > gene_start:
                    NovelExons_Subcategory.append(transcript + ",Beyond_First")
                elif novel_exon_start < gene_start and novel_exon_end > gene_end:
                    NovelExons_Subcategory.append(transcript + ",Internal_NovelExon")
                elif novel_exon_start < gene_end:
                    NovelExons_Subcategory.append(transcript + ",Beyond_Last")    
                elif novel_exon == "1":
                    NovelExons_Subcategory.append(transcript + ",First")
                else:
                    # check if the novel exon is the last exon (
                    # i.e. missing exon is beyond the final known exon in the transcript)
                    trans_pd = class_by_transcript_pd(transcript, All_FilteredParsed)
                    maxexondetected = max([int(i.replace("Exon","")) for i in trans_pd["TranscriptExon"].values])
                    if int(novel_exon) > maxexondetected:
                        NovelExons_Subcategory.append(transcript + ",Last")
                    else:
                        print("Transcript not classified for Novel Exons:", transcript)

        # Generate output table, of the number of exons that are beyond first, internal, beyond last for each transcript
        output_df = generate_split_table(NovelExons_Subcategory,"classification")
        
        # Write output of novel exon coordinates for each transcript with classification
        final = pd.concat([NE, output_df["classification"]], axis=1)
        final.to_csv(args.gene_stats_dir + args.genename + "_NE_coordinates.csv", index=False)


        # Generate list of transcripts with those classifications:
        # Prioritise Beyond_First Last, then Beyond First, Beyond  Last and Internal
        for transcript in output_df["transcriptID"].unique():
            novel_exon_types = list(output_df.loc[output_df["transcriptID"] == transcript,"classification"])
            if "Beyond_First_Last" in novel_exon_types:
                NExons_BeyondFirstLast.append(transcript)
            elif "Beyond_First" in novel_exon_types:
                NExons_BeyondFirst.append(transcript)
            elif "Beyond_Last" in novel_exon_types: 
                NExons_BeyondLast.append(transcript)
            elif "Internal_NovelExon":
                NExons_Internal.append(transcript)
            else:
                print("ERROR")

        # QC
        len(NExons_First) + len(NExons_Last) + len(NExons_Internal) + len(NExons_BeyondFirst) + len(NExons_BeyondLast) + len(NExons_BeyondFirstLast) == len(output_df["transcriptID"].unique())

    return output_df, NExons_BeyondFirst, NExons_BeyondFirstLast, NExons_BeyondLast, NExons_Internal, NExons_First, NExons_Last
        
