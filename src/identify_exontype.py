def identify_constitutive_alternative_exon(gene):
    
    '''
    Aim: Define the exons as either constitutive (in majority/all of the isoforms) or alternative
    Constitutive here is defined as exons that are in the majority of the isoforms 
    unless also specified in the UCSC genome browser basic annotation (i.e. some of the ANK1 shorter transcripts)
    1/ List the shorter transcripts not to be included for defining the gene exon classification
    otherwise majority of exons would be considered alternative even if present in all of the longer transcripts
    2/ Read in the flattened gene structure
    3/ Loop through each row, start and end coordinates and if not present within each exons of a set of transcripts 
        ==> alternative exon
    Prequisite: Require the gencode structure to be flattened and saved as gene_Manual_gencode.csv
    '''
    
    # Shorter transcripts not to be included in gene exon classification 
    # Mono-exons can complicate the classification and add exons that are unique, therefore specified for removal
    excluded_transcript = pd.read_csv("C:/Users/sl693/Desktop/Gencode_transcript_exclusion.csv")["excluded_transcript"].values
    monoexons = ["ENSMUST00000207525.1","ENSMUST00000207106.1","ENSMUST00000207596.1",'ENSMUST00000190478.1'] # Apoe, Picalm
    
    # Read in the gencode classification (previously manually curated to flatten the gene structure)
    gencode = pd.read_csv("C:/Users/sl693/Desktop/Ref/" + gene + "_Manual_gencode.csv") 
    gencode["start_end"] = gencode["start"].astype(str) + "_" + gencode["end"].astype(str)
    #print(gencode["transcript"].unique())
    
    # Remove monoexons to not loop through the unique exons for complication
    gencode = gencode.loc[~gencode["transcript"].isin(monoexons),]
    
    # Loop through each set of start and end coordinates
    # Then loop through each set of transcripts start and end coordinates
    alternative_exon = []
    for e in gencode["start_end"]:
        e_start = e.split("_")[0]
        e_end = e.split("_")[1]
        
        # loop through each transcript 
        for t in gencode["transcript"].unique(): 
            if t not in excluded_transcript: # not to include the short transcripts 
                trans_exon = gencode.loc[gencode["transcript"] == t,"start_end"].values
                trans_exon_start = [i.split("_")[0] for i in trans_exon]
                trans_exon_end = [i.split("_")[1] for i in trans_exon]
                # if the set of exon coordinates not within the transcript, or the start and end is not the same
                if e not in trans_exon and e_start not in trans_exon_start and e_end not in trans_exon_end:
                    alternative_exon.append(e)
    
    # obtain the unique updated exon number of the alternative exon
    updated_alt_exon = np.unique(gencode.loc[gencode["start_end"].isin(np.unique(alternative_exon))]["updated_exon_number"])
    
    # create final table 
    final = []
    for i in range(1, max(gencode["updated_exon_number"]) + 1):
        if i in updated_alt_exon:
            final.append(["Gencode_" + str(i),"alternative", gene])
        else:
            final.append(["Gencode_" + str(i),"constitutive", gene])
    
    return(final)