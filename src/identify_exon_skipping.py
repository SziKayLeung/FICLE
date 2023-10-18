#!/usr/bin/env python
# exon skipping

import pandas as pd
import numpy as np 
from collections import Counter
import sys

from prepare_and_parse import class_by_transcript_pd
from prepare_and_parse import class_by_transcript
from prepare_and_parse import generate_split_table
from prepare_and_parse import determine_order

def tabulate_exon_presence(gencode, df, All_FilteredParsed):
    df_transcript_id = df['transcript_id'].unique()

    '''
    Aim: On the parsed output i.e. for each transcript,
    1/ Iterate through each known gencode exon to identify if there is a transcript exon classified to it i.e not "No"
    2/ if there is a classification, ensure that it matches any of the labels except for intron retention
    :gencode = read gencode gtf

    Output: Table of all transcripts as rows and gencode exons as columns 
    1 = Gencode exon is present 
    0 = Gencode exon not present
    2 = Genocode exon is present as IR 
    1001 = Gencode exon is not registered*

    *Transcript is matching a reference transcript with fewer exons than other reference transcripts
    No need to address order (sense/antisense) as coordinates and exon naming 
    for gencode is addressed in determine_order() and parse_gencode_reference()
    for transcript df in parse_transcript()
    '''

    max_exon = max([int(i) for i in gencode["updated_exon_number"]])

    # All labels except IR
    Labels = ["AF", "MisMatch", "MisjumpMatch", "Match", "ExtendedA5","ExtendedA3", "TruncatedA5","TruncatedBothA3A5",
              "TruncatedA3","ExtendedBothA3A5"]
    Gencode_exons = ["Gencode_" + str(i) for i in range(1,max_exon+1)] 

    # Output as table 
    df2 = pd.DataFrame(columns = Gencode_exons)
    data = []

    # Iterate through each transcript
    for transcript in df_transcript_id:     

        tdat = class_by_transcript_pd(transcript, All_FilteredParsed)
        maxgencodexon = max(tdat["GencodeExon"])

        # Iterate through each gencode exon
        output = []
        for g in range(1,max_exon+1):
            # filter the classifcation for each corresponding gencode exon
            Gencode_filter = tdat.loc[tdat['GencodeExon'] == g,'Class'].values

            if g > max_exon:
                output.append(1001)
            elif len(list(filter(lambda x: "No" not in x, Gencode_filter))) != 0:
                if "IR" in Gencode_filter or "IRMatch" in Gencode_filter:
                    output.append(2) ## IR
                elif any(elem in Labels for elem in Gencode_filter):
                    output.append(1) ## 1 = Present
                else:
                    print("ERROR:", transcript)
                    sys.exit()
            else:
                output.append(0) # only No classifications                 

        # Create a dictionary of the output to append to the table
        zipped = zip(list(df2), output)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)

    df2 = df2.append(data, True)
    df2.index = df_transcript_id    

    ## Stats 
    # Count of the number of transcripts with number of exons (i.e 100 transcripts with 5 exons, 41 with 4 exons)
    num_exons_pertrans = (df2 == 1).astype(int).sum(axis=1)
    print("Number of transcripts with number of exons")
    print(Counter(num_exons_pertrans))
    
    return(df2)


'''
Aim: Used the parsed exon table of each transcript to identify exon skipping
1/ Iterate through each row of the exon_tab
2/ Iterate through each column of that row to identify the sequence of 0 and 1 for 
identifiying first exon, exon skipping, intron retention and if exon is present
:gencode = read gencode gtf
:exon_tab = table of exon presence output from tabulate_exon_presence()

** Intron Retention ** 
# IR defined by "2" until a "1" is present 
Important to distinguish exons that are IR rather than absent and thus mistakenly referred as ES 

** Absent Exons **
Some reference transcripts have fewer exons than other reference transcripts, 
thus a transcript can have fewer exons than all but is actually matching and does not exhibit ES
Reference exons that are missing in the transcript but are not classified w
ith a class (not even "No") given that the transcript exon is not present is captured as "1001"

** Exon Skipping **
# The "1" in the first sequence is always the first detected internal exon 
(note also other first exons but would be classified as novel)
# if 0 ==> "NA"

# For all other exons, i.e. not first or last exon 
# The "1" refers to present i.e. not skipping 
# The "0" refers to absent i.e skipping        

Output: Table of all transcripts as rows and gencode exons as columns 
1 = Gencode exon is present 
0 = Gencode exon not present
'''
def identify_exon_skipping(gencode,exon_tab,All_FilteredParsed):
        
    # Create gencode list of possible exon numbers
    max_exon = max([int(i) for i in gencode["updated_exon_number"]])
    Gencode_exons = ["Gencode_" + str(i) for i in range(1,max_exon+1)] 
        
    # empty dataframe for Exon skiipping
    ES = pd.DataFrame(columns = Gencode_exons)
    data = []
    
    # iterate through each row of exon table
    for index, row in exon_tab.iterrows():
        pre_value = []
        output = []
        
        # last exon detected for transcript
        dat = class_by_transcript_pd(index, All_FilteredParsed)
        transcript_lastexon = max([int(i) for i in dat["GencodeExon"].values])
        #print(transcript_lastexon)
      
        # iterate through each value in the sequence of that row
        # count = keeps score of the sequence value as iterating through row
        for count, ele in enumerate(row):

            # For the first exon, record it's first exon ("101") or alternative first exon ("NA")
            # if alternative first exon, any 0 after that does not refer to skipping but simply not present
            if count == 0:
                if ele == 1:
                    output.append("FirstExon")
                    pre_value = "101"
                elif ele == 2:
                    output.append("IR")
                    pre_value = "2"
                else:
                    output.append("NA")
                    pre_value = "NA" 
            
            # For internal exons
            elif count == transcript_lastexon - 1:
                if ele == 1:
                    output.append("No")
                else:
                    output.append("Yes")
            
            elif count > transcript_lastexon - 1:
                output.append("NA")
                
            elif 0 < count < max_exon - 1:
                if ele == 1001:
                    output.append("NA")
                    pre_value = "NA"
                if ele == 2:
                    output.append("IR")
                    pre_value = "2"
                elif ele == 1:
                    # if already detected first exon, then the 1 refers to present 
                    if pre_value == "101" or pre_value == "1" or pre_value == "0" or pre_value == "2": # not skipping
                        output.append("No")
                        pre_value = "1"
                    # if still not detected first exon, but 1 now, this refers to alterative first exon
                    elif pre_value == "NA":
                        output.append("FirstExon")
                        pre_value = "101"
                    else:
                        print("ERROR")
                else:
                    # if already detected first exon, then 0 refers to absent
                    if pre_value == "101" or pre_value == "1" or pre_value == "0": # Skipping
                        output.append("Yes") 
                        pre_value = "0"
                    elif pre_value == "2":
                        output.append("IR")
                        pre_value = "2"
                    elif pre_value == "NA": # Otherwise still NA, not present
                        output.append("NA")
                        pre_value = "NA"
                    else:
                        print("ERROR")

            else:
                if ele == 1:
                    output.append("Present_FinalExon")
                elif ele == 2:
                    output.append("IR_FinalExon")
                elif ele == 0:
                    output.append("Absent_FinalExon")
                else:
                    output.append("NA_FinalExon")

                

        zipped = zip(list(ES), output)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)        

    ES = ES.append(data, True)
    ES.index = exon_tab.index
    
    return(ES)

    
def gene_specific(gene, df, ES):
    
    def special_exon(special_start, special_end):
        captured_trans = []
        IRthreshold = 101
        # capture transcripts with first long exon, allowing wobble
        for index, row in df.iterrows():
            if abs(row["start"] - special_start) <= IRthreshold and abs(row["end"] - special_end) <= IRthreshold:
                captured_trans.append(row["transcript_id"])
        
        return(captured_trans)
        
    
    if gene == "Trem2":
        '''
        Situation: 1 reference transcript with long first exon that overlaps exon 1 and exon2 of other reference transcripts 
        Challenges: long-read transcripts mapped to this reference transcript appears be ES of Gencode 2
        &  undetected in categories with all exons matching due to only 4 rather than 5 exons
        Aim: Capture transcripts with this long exon, remodify ES tab as NA in Gencode 2, & seek all matching transcripts
        '''
        # long exon: chr17: 48346440 - 48348807
        long_1st = special_exon(48346440,48348807)
        print("Number of transcripts captured with the 1st exon:", len(long_1st))
        
        # Of these transcripts, isolate those with all matching exons to that of the reference 
        MatchLabels = ['Match_Gencode_1', 'Match_Gencode_3', 'Match_Gencode_4', 'Match_Gencode_5']
        AllKnownMatch_long1st = []
        for transcript in long_1st:
            class_transcript_exon,class_gencode_exon = class_by_transcript(transcript, All_FilteredParsed)
            if all(elem in class_gencode_exon for elem in MatchLabels):
                    AllKnownMatch_long1st.append(transcript)
        
        # Replace ES output Gencode 2 as NA
        ES.loc[long_1st, "Gencode_2"] = "NA"
        
        return AllKnownMatch_long1st, ES 

    if gene == "Cd33":
        '''
        Situation: 1 reference transcript with long first exon that overlaps exon 1 and exon2 of other reference transcripts 
        Challenges: long-read transcripts mapped to this reference transcript appears be ES of Gencode 2
        &  undetected in categories with all exons matching due to only 4 rather than 5 exons
        Aim: Capture transcripts with this long exon, remodify ES tab as NA in Gencode 2, & seek all matching transcripts
        '''
        # long exon: 43533058	43533101
        short_1st = special_exon(43533058,43533101)
        print("Number of transcripts captured with the 1st exon:", len(short_1st))
        
        # Of these transcripts, isolate those with all matching exons to that of the reference 
        MatchLabels = ['Match_Gencode_3', 'Match_Gencode_3', 'Match_Gencode_4', 'Match_Gencode_5','Match_Gencode_6',
                      'Match_Gencode_7','Match_Gencode_8']
        AllKnownMatch_short1st = []
        for transcript in short_1st:
            class_transcript_exon,class_gencode_exon = class_by_transcript(transcript, All_FilteredParsed)
            if all(elem in class_gencode_exon for elem in MatchLabels):
                    AllKnownMatch_short1st.append(transcript)
        
        # Replace ES output Gencode 2 as NA
        ES.loc[short_1st, "Gencode_1"] = "NA"
        ES.loc[short_1st, "Gencode_2"] = "NA"
        
        return AllKnownMatch_short1st, ES 
    
    if gene == "Clu":
        '''
        Situation: 3 reference transcripts with an upstream first exon and staggered exons 2, 3, etc 6 downstream
        Challenges: long-read transcripts with this upstream first exon but no exons 2-6 would be considered as exon skipping
        Aim: Capture transcripts with this alternative upstream first exons and exons 2-6, 
        remodify ES tab as NA in Gencode 2-6 if not classified as "FirstExon"
        '''
        manualgencode_output = manual_gencode()
        alt_1 = manualgencode_output[manualgencode_output["original_exon_number"] == 1]
        
        for index, row in alt_1.iterrows():
            Alt_first_exon = special_exon(row["start"], row["end"])
            for i in range(2,7): # Gencode 2 - 6
                for t in Alt_first_exon:
                    if ES.loc[t,"Gencode_" + str(i)] != "FirstExon":
                        ES.loc[t,"Gencode_" + str(i)] = "NA"       
        return(ES)
    
    if gene == "Apoe":
        manualgencode_output = manual_gencode()
        alt_1 = manualgencode_output[manualgencode_output["original_exon_number"] == 1]
        
        for index, row in alt_1.iterrows():
            Alt_first_exon = special_exon(row["start"], row["end"])
            for i in range(2,4): # Gencode 2 - 3
                for t in Alt_first_exon:
                    if ES.loc[t,"Gencode_" + str(i)] != "FirstExon":
                        ES.loc[t,"Gencode_" + str(i)] = "NA"
                
        return(ES)
    
    if gene == "Vgf":
        AllKnownMatch_short = []
        # Short transcript at Gencode Exon 4 therefore not AP if missing the first 3 exons
        for transcript in df["transcript_id"].unique(): 
            trans_df = class_by_transcript_pd(transcript, All_FilteredParsed)
            if len(trans_df[(trans_df["TranscriptExon"] == "Exon1") & (trans_df["GencodeExon"] == "Gencode_4")]) >= 1:
                AllKnownMatch_short.append(transcript)
        
        return AllKnownMatch_short
    
    if gene == "Snca":        
        Firstexon = []
        for index,row in ES.iterrows():
            if row["Gencode_1"] == "FirstExon" or row["Gencode_2"] == "FirstExon":
                Firstexon.append(index)
        ES.loc[Firstexon,"Gencode_6"] = "NA"
        return(ES)
    
    if gene == "Fyn":        
        Firstexon = []
        for index,row in ES.iterrows():
            if row["Gencode_1"] == "FirstExon" or row["Gencode_2"] == "FirstExon":
                Firstexon.append(index)
        ES.loc[Firstexon,"Gencode_4"] = "NA"
        ES.loc[Firstexon,"Gencode_5"] = "NA"
        ES.loc[Firstexon,"Gencode_6"] = "NA"
        return(ES)

    if gene == "Bin1":
        '''
        Situation: 1 reference transcripts with downstream first exon that is not present in ref transcripts with upstream first
        Challenges: long-read transcripts with upstream first exon but no exon 3 would be considered as exon skipping
        Aim: Capture transcripts with upstream first exons, remodify ES tab as NA in Gencode 3 if not classified as "FirstExon"
        '''
        Firstexon = special_exon(32377230, 32377490)
        
        for t in Firstexon:
            if ES.loc[t,"Gencode_3"] == "Yes": # if classified as Exon skipping
                ES.loc[t,"Gencode_3"] = "NA"                
        
        return(ES)

def exon_skip_not(ES_dat, lst, Gencode_exon):
        for t in lst:
            if ES_dat.loc[t,Gencode_exon] == "Yes": # if classified as Exon skipping
                ES_dat.loc[t,Gencode_exon] = "NA"

        return(ES_dat)
    
    
def skip_not_AFexons(ES, gencode):
    
    '''
    Aim: Do not include alternative first exons that are unique to only one transcript as exon skipping 
    Rationale: Inflate the number of ES events, which should reflect only exons that are constitutive or seen > 1
    1/ Use the gencode reference (modified and flattened) to extract unique exons (i.e. present only in one transcript)
    2/ If the unique exon refers to exon 1 from the transcript which is seen --> ES table modified to "NA" if "Yes"
    i.e not seen as ES
    
    An automated gene-specific way using the gencode structure flattened
    '''
    
    # find the maximum number of exons for looping downstream
    gencode = gencode.astype({"updated_exon_number": int})
    max_exon = max([i for i in gencode["updated_exon_number"]])
    assert max_exon > 0 
    
    # loop through each exon, count the number of times it is seen in the reference (num_entries)
    # if num_entries (number of times seen) is only once
    # check if that exon refers to exon 1 (original exon before flattening) 
    unique_AF_exon = []
    for num in range(1, max_exon):
        # tabulate the entries, which should be  more than 0 as flattened and seen
        num_entries = len(gencode[gencode["updated_exon_number"] == num])
        assert num_entries > 0
        # conditional that if one time
        if num_entries == 1:
            # and if it is exon 1
            if gencode.loc[gencode["updated_exon_number"] == num]["exon_number"].values[0] == "1":
                print("Alternative First and not Exon skipped: Exon", num)
                unique_AF_exon.append(num)
    
    # loop through the unique AF exons, and modify the ES table
    for i in unique_AF_exon:
        ES = exon_skip_not(ES, ES.index, "Gencode_" + str(i))
        
    return(ES)


def gene_specific_human(gene, df, ES, All_FilteredParsed, gencode_gtf):
    
    def special_exon(special_start, special_end):
        captured_trans = []
        IRthreshold = 101
        # capture transcripts with first long exon, allowing wobble
        for index, row in df.iterrows():
            if abs(row["start"] - special_start) <= IRthreshold and abs(row["end"] - special_end) <= IRthreshold:
                captured_trans.append(row["transcript_id"])
        
        return(captured_trans)
    
    def exon_skip_not(ES_dat, lst, Gencode_exon):
        for t in lst:
            if ES_dat.loc[t,Gencode_exon] == "Yes": # if classified as Exon skipping
                ES_dat.loc[t,Gencode_exon] = "NA"

        return(ES_dat)

    
    # human Gencode 3 is only present in minor isoform
    if gene == "ABCA1":
        ES = exon_skip_not(ES, ES.index, "Gencode_2")
        return(ES)
    
    elif gene == "ABCA7":
        ES = exon_skip_not(ES, ES.index, "Gencode_13")
        ES = exon_skip_not(ES, ES.index, "Gencode_47")
        return(ES)
    
    elif gene == "RHBDF2":
        ES = exon_skip_not(ES, ES.index, "Gencode_3")
        return(ES)
    
    elif gene == "MAPT":
        ES = exon_skip_not(ES, ES.index, "Gencode_2")
        ES = exon_skip_not(ES, ES.index, "Gencode_9")
        return(ES)
    
    elif gene == "CD33":
        ES = exon_skip_not(ES, ES.index, "Gencode_7")
        ES = exon_skip_not(ES, ES.index, "Gencode_8")
        return(ES)
    
    elif gene == "PTK2B":
        ES = exon_skip_not(ES, ES.index, "Gencode_8")
        ES = exon_skip_not(ES, ES.index, "Gencode_9")
        ES = exon_skip_not(ES, ES.index, "Gencode_10")
        return(ES)
    
    elif gene == "SORL1":
        ES = exon_skip_not(ES, ES.index, "Gencode_15")
        ES = exon_skip_not(ES, ES.index, "Gencode_23")
        ES = exon_skip_not(ES, ES.index, "Gencode_27")
        ES = exon_skip_not(ES, ES.index, "Gencode_34")
        return(ES)
    
    elif gene == "SNCA":
        ES = exon_skip_not(ES, ES.index, "Gencode_2")
        ES = exon_skip_not(ES, ES.index, "Gencode_3")
        ES = exon_skip_not(ES, ES.index, "Gencode_4")
        ES = exon_skip_not(ES, ES.index, "Gencode_5")
        return(ES)
    
    elif gene == "FYN":
        ES = exon_skip_not(ES, ES.index, "Gencode_2")
        ES = exon_skip_not(ES, ES.index, "Gencode_3")
        ES = exon_skip_not(ES, ES.index, "Gencode_5")
        ES = exon_skip_not(ES, ES.index, "Gencode_6")
        ES = exon_skip_not(ES, ES.index, "Gencode_7")
        ES = exon_skip_not(ES, ES.index, "Gencode_7")
        ES = exon_skip_not(ES, ES.index, "Gencode_9")
        ES = exon_skip_not(ES, ES.index, "Gencode_10")
        ES = exon_skip_not(ES, ES.index, "Gencode_11")
        ES = exon_skip_not(ES, ES.index, "Gencode_12")
        ES = exon_skip_not(ES, ES.index, "Gencode_13")
        return(ES)
    
    elif gene == "APP":
        ES = exon_skip_not(ES, ES.index, "Gencode_3")
        ES = exon_skip_not(ES, ES.index, "Gencode_4")
        return(ES)
    
    elif gene == "PICALM":
        ES = exon_skip_not(ES, ES.index, "Gencode_2")
        ES = exon_skip_not(ES, ES.index, "Gencode_3")
        ES = exon_skip_not(ES, ES.index, "Gencode_4")
        return(ES)

    elif gene == "CLU":
        ES = exon_skip_not(ES, ES.index, "Gencode_4")
        ES = exon_skip_not(ES, ES.index, "Gencode_5")
        ES = exon_skip_not(ES, ES.index, "Gencode_6")
        ES = exon_skip_not(ES, ES.index, "Gencode_7")
        return(ES)
    
    elif gene == "BIN1":
        ES = exon_skip_not(ES, ES.index, "Gencode_8")
        return(ES)
    
    else:
        print("Not performing gene specific annotations")
        return(ES)
    
    
def output_exon_skipping_stats(ES):
    
    '''
    Aim: Tabulate Stats for exon skipping 
    
    Output:
    ES_Count = Number of exons skipped per transcript
    ES_SpecificExonSkipped = Specific Gencode exon skipped per transcript
    ES_Transcripts = List of Transcripts with Exon Skipping 
    '''
    
    ES_Count = pd.DataFrame()
    ES_SpecificExonSkipped = pd.DataFrame()
    ES_Transcripts = []

    try:        
        # Table 1: Number of exons skipped per transcript
        # Count the number of "Yes" per row (Transcript) in ES table
        # Filter the rows (Transcript) with 0 exons skipped
        ES_Count = pd.DataFrame((ES.applymap(lambda x: str.count(x, 'Yes')) == 1).astype(int).sum(axis=1))
        ES_Count.columns = ['Count']
        ES_Count = ES_Count.loc[ES_Count["Count"] != 0,]

        # Table 2: Gencode Exon skipped per transcript
        output = []
        for transcript in ES_Count.index:
            # Filter the ES table by the transcript
            # Identify the column names with "Yes" and append to list
            ES_filtered = pd.DataFrame(ES.loc[transcript])
            ES_filtered_yes = ES_filtered[ES_filtered.iloc[:,0]=="Yes"].index.tolist()
            output.append([transcript + "," + s for s in ES_filtered_yes] )

        # Aggregate all list to one big list
        output_lst = sum(output, [])
        ES_SpecificExonSkipped = generate_split_table(output_lst,"ES")

        # Output 3: Total number of transcripts with that exon skipped
        # For each cell in ES, count the number of "Yes" i.e if gencode exon skipped for that transcript ==> 1
        # Count down the number of "1" for each column 
        num_exons_skipped = ES.applymap(lambda x: str.count(x, 'Yes'))
        num_exons_skipped_pertrans = (num_exons_skipped == 1).astype(int).sum(axis=0)
        #print("Number of transcripts with the gencode exon skipped, ignore first and last exon")
        #print(num_exons_skipped_pertrans)

        # Output 4: List of Transcripts with Exon Skipping 
        ES_Transcripts = list(ES_Count.loc[ES_Count["Count"] != 0,].index)

    except:
        print("No transcripts with exon skipping")
  
    
    return ES_Count, ES_SpecificExonSkipped, ES_Transcripts
    
    
    
    
