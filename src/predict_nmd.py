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


def split_ORF_id(t):
        # Split the ORF id in the blat analysis as includes the transcript id followed by the ORF reference 
        # Important for capturing the target id for downstream analysis 
        if t.count("_") == 2:
            s = t.split("_",2)[0]
        else:
            s = t.split("_",4)[0] + "_" + t.split("_",4)[1]
        return(s)
    
    
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