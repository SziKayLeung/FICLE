#!/usr/bin/env python
# coding: utf-8

# ### Table of Contents
# 
# * [Prepare reference and transcript gtf](#1)
# * [Matching Exons](#2)
# * [Tabulate Exons](#2b)
# * [Gene specific](#gene_specific)
# * [Exon Skipping](#3)
# * [Alternative Promoter](#4)
# * [Alternative Termination](#4b)
# * [Novel Exons](#5)
# * [Intron Retention](#6)
# * [Alternative A5 A3](#7)
# * [Alternative First](#8)
# * [Final Output](#9)
# * [Run](#Run)

# In[2]:


# Szi Kay Leung (sl693@exeter.ac.uk)


def identify_intron_retention_old(df, All_FilteredParsed):
    
    '''
    Aim: Tabulate the transcripts with intron retention from the All_FilteredParsed output
    For the counts, apply the function to disentangle the actual number of IR events from the number of exons with IR 
    1/ For transcripts with IR, order the number of exons with IR 
    i.e gencode_IR = [1, 2, 3, 5] ==> Exons 1, 3 and 5 classified with IR
    
    Output: 
    output_df: Table of Transcripts with corresponding exon with intron retention
    output_counts: The total number of transcripts with inton retention
    '''
    
    def demystify_IR(row):
        '''
        Demystify the number of IR events as the number of exons with IR != IR events
        :row = list of the gencode exon number classified with IR 
        if row = [1, 2, 3, 5] --> Exon 1, 2, 3 with IR = 3 exons with IR but only 1 IR event 
        if row = [1, 3, 4] --> same IR through Exons 3, 4 therefore total 2 IR events 
        
        # Catalogue the number of IR events based on the series of exons
        IR event is considered if the next exon is more than one up (i.e. not continuous)
        '''
        first = row[0] # first ordered gencode exon with IR
        prev = row[0]  # tabulate the previous exon
        count = 1 # count of IR, starting with 1
        for i in row:
            # if the exon number is not the first exon
            if(i != first):
                # if the exon number is more than the first exon + 1 i.e. a jump
                if(i != prev + 1):
                    # the include the count, and move the series along
                    count+=1
                # otherwise update the tally of previous count
                prev = i
        
        return(count)
    
    
    IR = []    
    IR_Count = []
    # Iterate through each transcript
    for transcript in df['transcript_id'].unique():
        class_transcript_exon,class_gencode_exon = class_by_transcript(transcript,  All_FilteredParsed)
        
        gencode_IR = []  # list of the gencode exon number classified with intron retention
        for i in class_gencode_exon:
            if "IR" in i:
                if i == "IR_Gencode_1": # do not include IR_Genocde --> treat that as A5A3
                    pass
                else:
                    IR.append(transcript + "," + i)
                    gencode_IR.append(int(i.split("_")[2])) # append the exon number
        
        if len(gencode_IR) != 0: # if there are exons with intron retention
            # change to integer, sort the order and apply the function
            gencode_IR = [int(item) for item in gencode_IR]
            gencode_IR = list(set(gencode_IR))
            gencode_IR.sort()
            IR_Count.append(transcript + "," + str(demystify_IR(gencode_IR)))
    
    
    output_df = pd.DataFrame()
    output_counts = pd.DataFrame()
    output_list = []
    try:
        output_df = generate_split_table(IR,"IR")
        output_counts = generate_split_table(IR_Count,"IR")
        output_list = output_df["transcript_id"].unique()
    except:
        print("No transcripts with intron retention")
    
    return output_df, output_counts, output_list


# In[26]:


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


# In[27]:


def identify_intron_retention(df, All_FilteredParsed, gencode):
    
    '''
    Aim: Tabulate the transcripts with intron retention from the All_FilteredParsed output
    For the counts, apply the function to disentangle the actual number of IR events from the number of exons with IR 
    1/ For transcripts with IR, order the number of exons with IR 
    i.e gencode_IR = [1, 2, 3, 5] ==> Exons 1, 3 and 5 classified with IR
    
    Last exon: if the IR is from the last exon of the updated reference
    
    Output: 
    output_df: Table of Transcripts with corresponding exon with intron retention
    output_counts: The total number of transcripts with inton retention
    '''
    
    def demystify_IR(row):
        '''
        Demystify the number of IR events as the number of exons with IR != IR events
        :row = list of the gencode exon number classified with IR 
        if row = [1, 2, 3, 5] --> Exon 1, 2, 3 with IR = 3 exons with IR but only 1 IR event 
        if row = [1, 3, 4] --> same IR through Exons 3, 4 therefore total 2 IR events 
        
        # Catalogue the number of IR events based on the series of exons
        IR event is considered if the next exon is more than one up (i.e. not continuous)
        '''
        first = row[0] # first ordered gencode exon with IR
        prev = row[0]  # tabulate the previous exon
        count = 1 # count of IR, starting with 1
        for i in row:
            # if the exon number is not the first exon
            if(i != first):
                # if the exon number is more than the first exon + 1 i.e. a jump
                if(i != prev + 1):
                    # the include the count, and move the series along
                    count+=1
                # otherwise update the tally of previous count
                prev = i
        
        return(count)
    
    max_exon = max([int(i) for i in gencode["updated_exon_number"]])
    
    IR = []    
    IR_Count = []
    IR_Exon1 = []
    IR_LastExon = []
    IRMatch_Count = 0
    IROnly_Count = 0
    # Iterate through each transcript
    #for transcript in ["PB.8560.273_TALONT001396616"]:
    for transcript in df['transcript_id'].unique():
        trans_df = class_by_transcript_pd(transcript, All_FilteredParsed)
        
        try: 
            gencode_IR = []  # list of the gencode exon number classified with intron retention
            IR_df = trans_df.loc[trans_df["Class"].isin(["IR","IRMatch"]),]
            IR.extend(transcript + "," + row["ClassGencodeExon"] for index,row in IR_df.iterrows())
            gencode_IR.extend(int(i.split("_")[1]) for i in IR_df["GencodeExon"].values) # append the exon number

            # last exon regardless of classification 
            IR_df["TranscriptExon"] = [int(i.replace("Exon","")) for i in IR_df["TranscriptExon"]]
            last_exon = IR_df.loc[IR_df['TranscriptExon'].idxmax(),"GencodeExon"].replace("Gencode_","")
            #print(last_exon)

            if len(IR_df["TranscriptExon"].unique()) == 1 and IR_df["TranscriptExon"].values[0] == 1:
                IR_Exon1.append(transcript)
            elif len(IR_df["TranscriptExon"].unique()) == 1 and IR_df["GencodeExon"].values[0] == "Gencode_" + str(max_exon):
                IR_LastExon.append(transcript)
            elif len(IR_df["TranscriptExon"].unique()) == 1 and IR_df["GencodeExon"].values[0] == "Gencode_" + str(last_exon): 
                IR_LastExon.append(transcript)
            else:
                pass
                   
            # to avoid counting Exon 1 and final exon as IR 
            IR = [i for i in IR if i.split(",",2)[1] not in ["IR_Gencode_1","IR_Gencode_" + str(max_exon),"IR_Gencode_" + str(last_exon)]]
            
            # to avoid double counting IRMatch         
            if "IRMatch" in IR_df["Class"].values: 
                IRMatch_Count = len(np.unique(IR_df[IR_df["Class"] == "IRMatch"]["TranscriptExon"].values))
                IRMatch_Gencode = [int(i.replace("Gencode_","")) for i in IR_df[IR_df["Class"] == "IRMatch"]["GencodeExon"].values]

                gencode_IR = [i for i in gencode_IR if i not in IRMatch_Gencode]

            if len(gencode_IR) != 0: # if there are exons with intron retention
                # change to integer, sort the order and apply the function
                gencode_IR = [int(item) for item in gencode_IR if item != 1]
                gencode_IR = list(set(gencode_IR))
                gencode_IR.sort()
                IROnly_Count = demystify_IR(gencode_IR)

            IR_Count.append(transcript + "," + str(IRMatch_Count + IROnly_Count))
        except:
            pass
    
    output_df = pd.DataFrame()
    output_counts = pd.DataFrame()
    output_list = []
    try:
        output_df = generate_split_table(IR,"IR")
        output_counts = generate_split_table(IR_Count,"IR")
        output_list = output_df["transcript_id"].unique()
    except:
        print("No transcripts with intron retention")
    
    return output_df, output_counts, output_list, IR_Exon1, IR_LastExon


# ### Alternative A5 A3 <a class="anchor" id="7"></a>

# In[28]:


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


# ### Alternative First <a class="anchor" id=8></a>

# In[29]:


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


# ### Final Output <a class="anchor" id="9"></a>

# In[30]:


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
    Transcript_Classifications.columns = ["Matching","MisandJumpMatch","MisMatch","MisjumpMatch", 
                                          "AllExons_A5A3","AllExons_MisMatch_A5A3","SomeMatch","MisjumpMatch_NotAll",
                                          "A5A3","AF","AP","AT","ES","IR","IR_Exon1Only", "IR_LastExonOnly",
                                          "NE_1st","NE_Int","NE_Last","NE_FirstLast"]
    Transcript_Classifications.index = df['transcript_id'].unique()
    
    # sum of novel exons across all novel exon cateogories
    Transcript_Classifications["NE_All"] = Transcript_Classifications[["NE_1st","NE_Int","NE_Last","NE_FirstLast"]].sum(axis=1)
    
    return(Transcript_Classifications)


# In[31]:


def prioritise_write_output(df, Transcript_Classifications, output_dir, gene, input_bed):
    
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
        elif row["MisandJumpMatch"] == 1:
            Output.append(transcript + "," + "MisandJumpMatch")
        elif row["MisjumpMatch"] == 1:
            Output.append(transcript + "," + "MisjumpMatch")
        elif row["MisMatch"] == 1:
            Output.append(transcript + "," + "MisMatch")
        elif row["Matching"] == 1:
            Output.append(transcript + "," + "Matching")
        elif row["AllExons_MisMatch_A5A3"] == 1:
            Output.append(transcript + "," + "AllExons_MisMatch_A5A3")
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
            Output.append(transcript + "," + "IR_ES")
        elif row["ES"] == 1 and row["NE_Int"] == 1 and row["AP"] == 0 and row["AT"] == 0:
            Output.append(transcript + "," + "ES_NeInt")
        elif row["NE_Int"] ==1 and row["AP"] == 0 or row["NE_Int"] ==1 and row["AT"] == 0: 
            # prioritise internal novel exons rather than novel exons from alternative promoter or termination
            Output.append(transcript + "," + "NE_Int")
        elif row["IR"] == 1:
            Output.append(transcript + "," + "IR")
        elif row["ES"] == 1:
            Output.append(transcript + "," + "ES")
        elif row["SomeMatch"] == 1 and row["AP"] == 0 or row["SomeMatch"] == 1 and row["AT"] == 0:
            Output.append(transcript + "," + "SomeMatch")
        elif row["MisjumpMatch_NotAll"] == 1:
            Output.append(transcript + "," + "MisjumpMatch_NotAll")
        elif row["AP"] == 1:
            Output.append(transcript + "," + "AP")
        elif row["AT"] == 1:
            Output.append(transcript + "," + "AT")
        elif row["AF"] == 1:
            Output.append(transcript + "," + "AF")
        elif row["AllExons_A5A3"] == 1:
            Output.append(transcript + "," + "AllExons_A5A3")
        elif row["A5A3"] == 1:
            Output.append(transcript + "," + "A5A3")
        else:
            print("Not Classified for Final output:" + transcript)
    
    # split coloured sorted bed file 
    bed = pd.read_csv(input_bed, sep = "\t", header = None)
    
    # write the output to a log file 
    output_file = open(output_dir + "/" + gene + "_" + "Locator_Bedfiles.txt","w")
    for element in Output:
        output_file.write(element + "\n")  

    
    # Create tuple for quick access for generating output files
    Final_Output = [(i.split(",",2)[0], i.split(",",2)[1]) for i in Output]
    for i in set([x[1] for x in Final_Output]):
        path = output_dir + "/" + gene + "_" +  str(i) 
        print("Filtering for:" + i)
        filter_output = [x[0] for x in list(filter(lambda cate: cate[1] == i, Final_Output))]
        # subset the bedfiles from the group of transcript per categories
        group = bed[bed[3].isin(filter_output)]
        print(len(group.index))
        # Splitting into multiple files by 1000 if more than 100 transcripts
        if len(group) > 1000:
            for count, i in enumerate(np.array_split(group, np.ceil(len(group)/1000)),1):
                print("Splitting and writing out to file", count)
                i.to_csv(path + "_" + str(count) + "_sorted_coloured.bed12", sep = "\t", index = None, header = None)
        else:
            group.to_csv(path + "_sorted_coloured.bed12", sep = "\t", index = None, header = None)
   


# In[32]:


def populate_classification(Transcript_Classifications, A5A3, IR_Counts, ES_Count, NE, NE_pertrans_classify_counts):
    '''
    Aim: Repopulate the "1" in the classification table with the actual number of events per transcript
    Using the counts stats output across the event types 
    1/ A5A3 = Number of exons with alternative 5' start or alternative 3' end sites - defined by extended or truncated 
    2/ ES = Number of exons skipped 
    3/ NE... = Number of novel exons
    '''
        
    def generate_NE_dict(cate):
        dat = NE_pertrans_classify_counts[NE_pertrans_classify_counts["novelexons"] == cate]
        dat_dict = dict(zip(dat["transcript_id"],dat[0]))
        return(dat_dict)
    
    def remap_transcript_classification(col, input_dict):
        Transcript_Classifications['isoform'] = Transcript_Classifications.index
        new_col = Transcript_Classifications['isoform'].map(input_dict).fillna(Transcript_Classifications[col])
        return(new_col)
    
    
    if(sum(Transcript_Classifications["A5A3"]) > 0): Transcript_Classifications[['A5A3']] = remap_transcript_classification("A5A3", collections.Counter(A5A3["transcript_id"]))
    if(sum(Transcript_Classifications["IR"]) > 0): Transcript_Classifications[['IR']] = remap_transcript_classification("IR", dict(zip(IR_Counts["transcript_id"],IR_Counts["IR"])))
    Transcript_Classifications['ES'] = remap_transcript_classification("ES", dict(zip(ES_Count.index,ES_Count["Count"])))
    
    if(len(NE) > 0):
        Transcript_Classifications['NE_1st'] = remap_transcript_classification("NE_1st", generate_NE_dict("Beyond_First"))
        Transcript_Classifications['NE_Int'] = remap_transcript_classification("NE_Int", generate_NE_dict("Internal_NovelExon"))
        Transcript_Classifications['NE_Last'] = remap_transcript_classification("NE_Last", generate_NE_dict("Beyond_Last"))
        Transcript_Classifications['NE_FirstLast'] = remap_transcript_classification("NE_Last", generate_NE_dict("Beyond_First_Last"))
    
    return(Transcript_Classifications)


# In[33]:


def delete_make_dir(output_dir):
    try:
        shutil.rmtree(output_dir)
    except:
        print("1st Round of annotations")
    os.mkdir(output_dir)
    os.mkdir(output_dir + "/Stats")


# ### NMD <a class="anchor" id="10"></a>

# In[34]:


def predict_nmd(gencode, transcript):
    
    '''
    Aim: Predict nonsense-mediated decay for transcript using: 
    - ORF prediction from CPAT --> ORF length, start and end
    - Blat analysis of Best ORF.fa from CPAT to mouse reference genome --> genomic coordinates of CDS
    
    ## Assuming if there is only one entry for the blat analysis for the transcript and in the correct chromosome:
    1/ Extract basic information about the ORF from the CPAT output
    2/ Extract the genomic coordinates of the ORF from the Blat analysis (start and end) 
    start = A from AUG and end is the last position including the STOP codon
    3/ Extract the coordinte of the last exon-exon junction = start coordinate of the last exon
    4/ Determine the distance between the ORF end genomic coordinate and the last junction
    5/ Determine coding status from CPAT output (mouse threshold: 0.44)
    6/ Deduce NMD status based on the distance and if coding 
    
    # depends on gene orientation
    sense: 
              |----------->--------| Last Exon ---------------| 
        ORF---|                                                    Predicted_NMD = Yes 
        ORF---|----------->----------------|                       Predicted_NMD = No
    distance between start coordinate of last exon and end genomic coordinate of ORF
    
    antisense:
       |-----------Last Exon | ----------<---------|
                                                   |-----ORF       Predicted_NMD = Yes
                 |-----------| ----------<---------|-----ORF       Predicted_NMD = No
    distance between end coordinate of the first exon (due to antisense) and start genomic coordinate of ORF    
    '''
    
    # chromosome of the target gene for later validation
    chrom = df.loc[df["transcript_id"] == transcript,"seqname"].values[0]
    
    def split_ORF_id(t):
        # Split the ORF id in the blat analysis as includes the transcript id followed by the ORF reference 
        # Important for capturing the target id for downstream analysis 
        if t.count("_") == 2:
            s = t.split("_",2)[0]
        else:
            s = t.split("_",4)[0] + "_" + t.split("_",4)[1]
        return(s)
    
    # split the transcript id in the blat analysis and extract the transcript of interest 
    blast_orf["transcript_id_split"] = [split_ORF_id(i) for i in blast_orf["transcript_id"]]
    tb_orf = blast_orf.loc[blast_orf["transcript_id_split"] == transcript,]
    
    # sense or antisense
    gencode, order = determine_order(gencode)

    # if there is only one output from the blat analysis corresponding to the transcript and in the correct chromosome
    if(len(tb_orf) == 0):
        return([transcript, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL"])
    else:
        if(len(tb_orf) >= 1):
            #print("More than one blat hit:", transcript)
            # selecting the blat hit with the greatest match
            tb_orf = tb_orf.loc[tb_orf["Tname"] == chrom,]
            tb_orf = tb_orf.loc[tb_orf['match'].idxmax(),]
        elif str(tb_orf["Tname"].values[0]) != chrom:
            print("Wrong chromosome for NMD:", transcript)
            sys.exit(0)
        else:
            pass

        # ORF start, end and length
        tORF = ORF.loc[ORF["seq_ID"] == transcript,]
        ORF_start = int(tORF["ORF_start"])
        ORF_end = int(tORF["ORF_end"])
        ORF_length = int(tORF["ORF"])

        # genomic position 
        cds_genomic_start = int(tb_orf["Tstart"])
        cds_genomic_end = int(tb_orf["Tend"])

        # last exon-exon junction
        # determine distance to last junction
        if order == "sense":
            lastjunc=df.loc[df["transcript_id"] == transcript,"start"].values[-1]
            dist_to_lastjunc = cds_genomic_end - lastjunc
        else: 
            lastjunc=df.loc[df["transcript_id"] == transcript,"end"].values[0]
            dist_to_lastjunc = lastjunc - cds_genomic_start

        # determine coding status
        coding_status = "Coding" if float(tORF["Coding_prob"]) >= 0.44 else "Non_Coding"

        # NMD status
        if coding_status == "Coding":
            if (dist_to_lastjunc) < -50:
                NMD_status = "Yes"
            else:
                NMD_status = "No"
        else:
            NMD_status = "NA"
    
    return([transcript, ORF_start, ORF_end, ORF_length, cds_genomic_start, cds_genomic_end, lastjunc, dist_to_lastjunc, coding_status, NMD_status])


# In[35]:


def split_ORF_id(t):
        # Split the ORF id in the blat analysis as includes the transcript id followed by the ORF reference 
        # Important for capturing the target id for downstream analysis 
        if t.count("_") == 2:
            s = t.split("_",2)[0]
        else:
            s = t.split("_",4)[0] + "_" + t.split("_",4)[1]
        return(s)


# In[36]:


def call_NMD_prediction():
    
    '''
    Aim: Loop through the NMD prediction function through all the transcripts associated with target gene 
    output: Large dataframe with all the transcripts and ORF info
    '''
    NMD_predict_output = []
    
    for transcript in np.unique(df["transcript_id"].values):
        NMD_predict_output.append(predict_nmd(gencode, transcript))
        
    Final = pd.DataFrame(NMD_predict_output)
    Final.columns = ["transcript","ORF_start","ORF_end","ORF_length","cds_genomic_start","cds_genomic_end","lastjunc","dist2lastjunc","Coding_status","NMD_status"]
    return(Final)


# ### Run <a class="anchor" id="Run"></a>

# In[37]:


def main(gene):
    root_dir = "C:/Users/sl693/Dropbox/UpdatedAnnotations/ADBDR/"
    input_bed = root_dir + "IsoSeqONT_final_genename_sorted_coloured.bed12"
    gencode_gtf = "C:/Users/sl693/Dropbox/Annotations/Ref/" + gene + "_gencode.gtf"
    input_gtf = "C:/Users/sl693/Dropbox/UpdatedAnnotations/IsoSeqONT_final_genename_corrected.gtf"
    noISM_path = "C:/Users/sl693/Dropbox/UpdatedAnnotations/IsoSeqONT_final_genename_classification_noISM.txt"
    output_dir = "C:/Users/sl693/Dropbox/UpdatedAnnotations/" + gene
    output_log = root_output_dir + gene + "_parsed_transcripts.txt" 
    output_log = open(output_log, "w")
    
    # ORF for NMD Prediction
    ORF = pd.read_csv(root_dir + "ADBDR_cpat.ORF_prob.best.tsv", sep= "\t")
    #blast_orf = pd.read_csv("C:/Users/sl693/Desktop/blat.outfinal.txt", sep= "\t")

    # prepare directory
    delete_make_dir(output_dir)

    ## read in gencode and transcriptome gtf
    try: 
        gencode = pd.read_csv("C:/Users/sl693/Desktop/Ref/" + gene + "_Manual_gencode.csv") 
        gencode.rename({'transcript': 'transcript_id'}, axis=1, inplace=True)
    except: 
        gencode = parse_gencode_reference()
    df = parse_transcriptome_gtf(input_gtf, gene, noISM_path)  
    
    # Parse through the transcriptome, classify and filter 
    All_FilteredParsed = []
    for count, transcript in enumerate(df['transcript_id'].unique()):
        #print(transcript)
        parsed = parse_transcript(gencode,df,transcript,10, 100)
        All_FilteredParsed.append(filter_parsed_transcript(gene,gencode,parsed, output_log))
        if count%50==0:
            print("Parsing through transcript", count)

    # Aggregate all the filtered parsed output from each transcript into one big list
    All_FilteredParsed = [x for l in All_FilteredParsed for x in l]

    # QC: Check that the transcripts in the original transcriptome gtf captured in the final big list
    if set(df["transcript_id"].unique()) != set(set([i.split(';',3)[0] for i in All_FilteredParsed])):
        sys.exit(-1)

    #All_Detected_Exons_Matching = identify_all_matching_exons(gencode, df, All_FilteredParsed)
    AllKnownMatch, Mis2Match, MisMatch, MisjumpMatch, OtherClass, OtherClassMis2Match, SomeMatch, MisjumpMatch_NotAll = identify_all_matching_exons(gencode, df, All_FilteredParsed)
    #All_Detected_Exons = identify_all_exons(gencode, df, All_FilteredParsed)
    #All_Detected_Exons_NonMatching = list(set(All_Detected_Exons) - set(All_Detected_Exons_Matching))

    # Tabulate exon presence 
    print("Tabulating exon presence")
    exon_tab = tabulate_exon_presence(gencode, df, All_FilteredParsed)

    # Exon Skipping 
    print("Processing transcripts for exon skipping")
    ES = identify_exon_skipping(gencode,exon_tab)
    if gene == "Trem2": 
        AllKnownMatch_long1st, ES = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)
        AllKnownMatch = AllKnownMatch + AllKnownMatch_long1st

    if gene == "Cd33": 
        AllKnownMatch_long1st, ES = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)
        AllKnownMatch = AllKnownMatch + AllKnownMatch_long1st

    if gene == "Clu": 
        ES = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)

    if gene == "Apoe": ES = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)
    if gene == "Snca": ES = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)
    if gene == "Fyn": ES = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)
    if gene == "Ank1": ES = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)
    if gene == "Sorl1": ES = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)
        

    ES_Count, ES_SpecificExonSkipped, ES_Transcripts = output_exon_skipping_stats(ES)

    ## QC
    #Gencode_exons = ["Gencode_" + str(i) for i in range(1,5+1)] 
    #ZeroExonSkipping = list(exon_tab[exon_tab[Gencode_exons].eq(1).all(1)].index)
    #set(ZeroExonSkipping) == set(All_Detected_Exons)

    # Alternative First Promoter     
    Alternative_First_Promoter = identify_alternative_promoter(df, ES,gene, gencode, gencode_gtf, All_FilteredParsed)
    if gene == "Vgf": 
        AllKnownMatch_short = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)
        Alternative_First_Promoter = list(set(Alternative_First_Promoter) - set(AllKnownMatch_short))

    # Alternative Termination 
    Alternative_Termination = identify_alternative_termination(df, gencode, All_FilteredParsed)

    # Alternative First Exon 
    Alternative_First_exon = identify_alternative_first(df, All_FilteredParsed)

    # Novel Exons 
    print("Identifying transcripts with novel exons")
    NE, NE_novel_co = identify_novel_exon(df, gencode, All_FilteredParsed)
    NE_classify, NExons_BeyondFirst, NExons_BeyondFirstLast, NExons_BeyondLast, NExons_Internal = classify_novel_exon(gencode, df, NE)
    NE_pertrans_classify_counts = pd.DataFrame()
    if len(NE) > 0: NE_pertrans_counts, NE_classify_counts, NE_pertrans_classify_counts = novel_exon_stats(NE, NE_classify)

    # Intron Retention 
    print("Identifying transcripts with intron retention")
    IR, IR_Counts, IR_Transcripts, IR_Exon1, IR_LastExon = identify_intron_retention(df, All_FilteredParsed, gencode)

    # Alternative A5' and A3' 
    print("Identifying transcripts with alternative 5' and 3' sites")
    A5A3, A5A3_pertrans_counts, A5A3_Counts, A5A3_Transcripts = identify_A5A3_transcripts(df, All_FilteredParsed)

    # Final Output 
    # Event lists = each category is generated from previous functions and contain list of transcripts under that event type
    categories = [AllKnownMatch, Mis2Match, MisMatch, MisjumpMatch, OtherClass, OtherClassMis2Match, SomeMatch,MisjumpMatch_NotAll,
                  A5A3_Transcripts, Alternative_First_exon, Alternative_First_Promoter, Alternative_Termination,
                  ES_Transcripts, IR_Transcripts, IR_Exon1, IR_LastExon, NExons_BeyondFirst, NExons_Internal,NExons_BeyondLast,NExons_BeyondFirstLast]
    
    Transcript_Classifications = generate_aggregated_classification(df, categories)
    prioritise_write_output(df, Transcript_Classifications, output_dir, gene, input_bed)
    Transcript_Classifications_Remapped = populate_classification(Transcript_Classifications, A5A3, IR_Counts, ES_Count, NE, NE_pertrans_classify_counts)


    # All other stats 
    exon_tab.to_csv(output_dir + "/Stats/" + gene + "_Exon_tab.csv")
    ES.to_csv(output_dir + "/Stats/" + gene + "_Exonskipping_generaltab.csv")
    A5A3.to_csv(output_dir + "/Stats/" + gene + "_A5A3_tab.csv")
    IR.to_csv(output_dir + "/Stats/" + gene + "_IntronRetention_tab.csv", index = False)
    ES_SpecificExonSkipped.to_csv(output_dir + "/Stats/" + gene + "_Exonskipping_tab.csv")
    Transcript_Classifications_Remapped.to_csv(output_dir + "/Stats/" + gene + "_Final_Transcript_Classifications.csv")
    gencode.to_csv(output_dir + "/Stats/" + gene + "flattened_gencode.csv")
    print("All Done!")

