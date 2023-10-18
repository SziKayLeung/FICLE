#!/usr/bin/env python
# finalise output

import numpy as np
import pandas as pd
import collections
import os


def generate_aggregated_classification(df, categories):
    
    '''
    Aim: Create one big table with a row for each transcript and a "1" if the event is present 
    Table is for downstream prioritising of the presence of event and for visualisation 
    1/ Create empty array and populate if the transcript is recorded within the event lists
    '''
    
    # Create empty matrix
    matrix = np.zeros((len(df['transcript_id'].unique()),len(categories)), dtype=int)
    
    # loop through the transcripts, and populate matrix
    for row, transcript in enumerate(df['transcript_id'].unique()):     
        for col, cat in enumerate(categories):
            # if transcript is in the event list, convert 0 to 1 in matrix
            if transcript in cat:
                matrix[row, col] = 1


    Transcript_Classifications = pd.DataFrame(matrix)
    Transcript_Classifications.columns = ["Matching","SomeMatch",
                                          "A5A3","AF","AP","AT","ES","IR","IR_Exon1Only", "IR_LastExonOnly",
                                          "NE_1st","NE_Int","NE_Last","NE_FirstLast"]
    Transcript_Classifications.index = df['transcript_id'].unique()
    # note just for prioritise_write_output as need to recreate column for exact number of novel exons 
    Transcript_Classifications["NE_All"] = Transcript_Classifications[["NE_1st","NE_Int","NE_Last","NE_FirstLast"]].sum(axis=1)
        
    return(Transcript_Classifications)



def prioritise_write_output(args, df, Transcript_Classifications):
    
    '''
    Aim: Generate gtf for the different categories, prioiritised by: 
    1) Matching, 2) AF, 3) novel exon first, 4) novel exon last, 5) novel exon first and last
    6) IR and ES and Novel exon internal 7) IR and ES, 8) ES and Novel exon internal 
    9) Novel exon internal, 10) IR, 11) ES, and 12) A5, A3
    
    Note: Multiple transcripts with multiple event types, however visaulisation based on prioritisation of event types
    A5A3 not included in the list as transcripts under this category will be captured in others
    
    :Transcript_classifications: output table generated from generate_aggregated_classification()
    :output_dir: path for output gtf
    :gene: target gene of interest for annotations (prefix for output gtf)
    '''
    
    Output = []

    for transcript in Transcript_Classifications.index:
        row = Transcript_Classifications.loc[transcript]
        if row["NE_All"] != 0 and row["Matching"] == 1:
            Output.append(transcript + "," + "NEMatching")
        elif row["Matching"] == 1:
            Output.append(transcript + "," + "Matching")
        elif row["NE_1st"] ==1:
            Output.append(transcript + "," + "NE_First")
        elif row["NE_Last"] == 1:
            Output.append(transcript + "," + "NE_Last")
        elif row["NE_FirstLast"] == 1:
            Output.append(transcript + "," + "NE_firstLast")
        elif row["IR_Exon1Only"] == 1 and row["AP"] == 1:
            Output.append(transcript + "," + "IRExon1_AP") 
        elif row["IR_Exon1Only"] == 1 and row["AP"] == 0:
            Output.append(transcript + "," + "IR_FirstExonOnly")
        elif row["IR_LastExonOnly"] == 1 and row["ES"] == 0:
            Output.append(transcript + "," + "IR_LastExonOnly")
        elif row["IR"] == 1 and row["ES"] == 1 and row["NE_Int"] == 1 and row["AP"] == 0 and row["AT"] == 0:
            Output.append(transcript + "," + "IR_ES_NEInt") 
        elif row["IR"] == 1 and row["ES"] == 1:
            Output.append(transcript + "," + "IR_ES_Both")
        elif row["ES"] == 1 and row["NE_Int"] == 1 and row["AP"] == 0 and row["AT"] == 0:
            Output.append(transcript + "," + "ES_NeInt_Both")
        elif row["NE_Int"] ==1 and row["AP"] == 0 or row["NE_Int"] ==1 and row["AT"] == 0: 
            # prioritise internal novel exons rather than novel exons from alternative promoter or termination
            Output.append(transcript + "," + "NEIntOnly")
        elif row["IR"] == 1:
            Output.append(transcript + "," + "IROnly")
        elif row["ES"] == 1:
            Output.append(transcript + "," + "ESOnly")
        elif row["SomeMatch"] == 1 and row["AP"] == 0 or row["SomeMatch"] == 1 and row["AT"] == 0:
            Output.append(transcript + "," + "SomeMatch")
        elif row["AP"] == 1:
            Output.append(transcript + "," + "APOnly")
        elif row["AT"] == 1:
            Output.append(transcript + "," + "AT")
        elif row["AF"] == 1:
            Output.append(transcript + "," + "AF")
        elif row["A5A3"] == 1:
            Output.append(transcript + "," + "A5A3")
        else:
            print("Not Classified for Final output:" + transcript)
    
    # split coloured sorted bed file 
    bed = pd.read_csv(args.input_bed, sep = "\t", header = None)
    
    # write the output to a log file 
    output_file = open(args.gene_tracks_dir + "/" + args.genename + "_" + "Locator_Bedfiles.txt","w")
    for element in Output:
        output_file.write(element + "\n")  

    
    # Create tuple for quick access for generating output files
    Final_Output = [(i.split(",",2)[0], i.split(",",2)[1]) for i in Output]
    print("**** Final classification of transcripts by...")
    for i in set([x[1] for x in Final_Output]):
        path = args.gene_tracks_dir + "/" + args.genename + "_" +  str(i) 
        filter_output = [x[0] for x in list(filter(lambda cate: cate[1] == i, Final_Output))]
        # subset the bedfiles from the group of transcript per categories
        group = bed[bed[3].isin(filter_output)]
        print(i,":", len(group.index))
        # Splitting into multiple files by 1000 if more than 100 transcripts
        if len(group) > 1000:
            for count, i in enumerate(np.array_split(group, np.ceil(len(group)/1000)),1):
                print("Splitting and writing out to file", count)
                i.to_csv(path + "_" + str(count) + "_sorted_coloured.bed12", sep = "\t", index = None, header = None)
        else:
            group.to_csv(path + "_sorted_coloured.bed12", sep = "\t", index = None, header = None)
   
   
def populate_classification(args, Transcript_Classifications, A5A3_Counts, IR_Counts, ES_Counts, NE_Counts):
    '''
    Aim: Repopulate the "1" in the classification table with the actual number of events per transcript
    Using the counts stats output across the event types 
    1/ A5A3 = Number of exons with alternative 5' start or alternative 3' end sites - defined by extended or truncated 
    2/ ES = Number of exons skipped 
    3/ NE... = Number of novel exons
    '''
        
    def generate_NE_dict(cate):
        dat = NE_Counts[NE_Counts["NEtype"] == cate]
        dat_dict = dict(zip(dat["transcriptId"],dat["numNovelExons"]))
        return(dat_dict)
    
    def remap_transcript_classification(col, input_dict):
        Transcript_Classifications['isoform'] = Transcript_Classifications.index
        new_col = Transcript_Classifications['isoform'].map(input_dict).fillna(Transcript_Classifications[col])
        return(new_col)
    
   
    if(sum(Transcript_Classifications["A5A3"]) > 0): Transcript_Classifications['A5A3'] = remap_transcript_classification("A5A3", collections.Counter(A5A3_Counts["transcriptID"]))
    if(sum(Transcript_Classifications["IR"]) > 0): Transcript_Classifications['IR'] = remap_transcript_classification("IR", dict(zip(IR_Counts["transcriptID"],IR_Counts["numEvents"])))
    Transcript_Classifications['ES'] = remap_transcript_classification("ES", dict(zip(ES_Counts.index,ES_Counts["numEvents"])))
     
    if(len(NE_Counts) > 0):
        Transcript_Classifications['NE_1st'] = remap_transcript_classification("NE_1st", generate_NE_dict("Beyond_First"))
        Transcript_Classifications['NE_Int'] = remap_transcript_classification("NE_Int", generate_NE_dict("Internal_NovelExon"))
        Transcript_Classifications['NE_Last'] = remap_transcript_classification("NE_Last", generate_NE_dict("Beyond_Last"))
        Transcript_Classifications['NE_FirstLast'] = remap_transcript_classification("NE_Last", generate_NE_dict("Beyond_First_Last"))
        # sum of novel exons after remapping across all novel exon cateogories
        Transcript_Classifications["NE_All"] = Transcript_Classifications[["NE_1st","NE_Int","NE_Last","NE_FirstLast"]].sum(axis=1)
        
    # write output 
    Transcript_Classifications.drop('isoform', axis=1, inplace=True)
    Transcript_Classifications.index.name = 'isoform'
    Transcript_Classifications = Transcript_Classifications.astype(int)
    Transcript_Classifications.to_csv(args.gene_stats_dir + args.genename + "_final_transcript_classifications.csv")
    
    return(Transcript_Classifications) 