import pandas as pd
import numpy as np
from prepare_and_parse import reidentify_isoform_dataset

def identify_mapt_exon2_3_17(species, gencode):
    '''
    Aim: Identify from the flattened structure of MAPT (using the parse_gencode_reference()), the exons corresponding to
    exon 2, exon 3 and exon 10 in the genome ("species")
    '''
    if species == "human":
        # <exon number> <exon start coordinate> <exon end coordinate>
        # coordinates from hg38, manually checked across literature
        exon2_co = ["2","45971859","45971945"]
        exon3_co = ["3","45974385","45974471"]
        exon10_co = ["10","46010310","46010402"]
    elif species == "mouse":
        exon2_co = ["2","104287123","104287209"]
        exon3_co = ["3","104289908","104289994"]
        exon10_co = ["10","104318157","104318249"]
   
    # function to match the exon coordinates against the reference parsed genome database to find the updated exon number
    # updated exon number refers to the exon number after flattening the gene structure
    def extract_mapt_equivalent_exon(exon_co):
        # extract the updated exon number
        mod_exon = np.unique(gencode.loc[(gencode["start"] == int(exon_co[1])) & 
                                         (gencode["end"] == int(exon_co[2])),"updated_exon_number"])
        # should only have one input 
        assert len(mod_exon) == 1
        # return exon number and the modified/updated exon number 
        return(int(exon_co[0]), mod_exon[0])
    
    # extract the modified exon number for exons 2, 3, and 10    
    output = []
    for i in [exon2_co, exon3_co, exon10_co]:
        output.append(extract_mapt_equivalent_exon(i))
    output = pd.DataFrame(output, columns = ["MAPT_orig_exon","MAPT_mod_exon"])
    
    # Add the word "Gencode" in front of the updated exon number for downstream subsetting
    output['MAPT_mod_exon'] = 'Gencode_' + output['MAPT_mod_exon'].astype(str)
    
    return(output)


def identify_first_exon_pos(lst):
    '''
    Aim: Find the position of the first exon i.e the position at which it's no longer "0" in the exon_tab
    Params:
        lst = list of 0s and 1s from the row of the exon table
    '''
    # loop through the list
    for count, exon in enumerate(lst):
        # continue through the list until it is no longer "1"
        if exon != 0:
            # + 1 as the count starts from 0 rather than 1, preserve exon position
            first_exon_occurence = count + 1
            break
        else:
            continue
    return(first_exon_occurence)


def classify_mapt_isoforms(species, mapt_exon_tab, gencode):
    '''
    Aim: Classify the mapt isoforms by the the skipping and presence of exons 2, 3, and 10 (according to literature)
    params: 
        species = human/mouse; this specifies which coordinates of exons 2, 3 and 10 to use in identify_mapt_exon2_3_17()
        mapt_exon_tab = table of the presence/absence of exons 
        
    The mapt nomenclature is defined by the presence/absence of exons 2, 3 and 10; 
    Exon 2 and 3 presence determines whether it's 0N, 1N or 2N (where they are both skipped)
    Exon 10 presence determines whether's its 3R or 4R (where it is included)
    '''
    
    # apply the function to identify the updated exon number to work with 
    mapt_exons = identify_mapt_exon2_3_17(species, gencode)
   
    # to identify transcripts that have an alternative first exon downstream of exons 2 and exon 3
    missing_e2_e3 = []
    for index, row in mapt_exon_tab.iterrows():
        first_exon_occurence = identify_first_exon_pos(row.values) 
        if first_exon_occurence > 3:
            missing_e2_e3.append(index)
    
    # extract the columns based on the updated exon number    
    df = mapt_exon_tab[mapt_exon_tab.columns.intersection(mapt_exons["MAPT_mod_exon"])]
    
    output = []
    # loop through each row and append the classification
    for index, row in df.iterrows():
        # if the transcript is within the list of missing exons 2 and exon 3, but exon 10 present
        if index in missing_e2_e3 and row[2] == 1:
            output.append("E2E3'4R")
        # if the transcript is within the list of missing exons 2 and exon 3, but exon 10 skipped
        elif index in missing_e2_e3 and row[2] == 0: 
            output.append("E2E3'3R")
        else:
            # exon 2, exon 3 skipped, exon 10 present
            if row[0] == 0 and row[1] == 0 and row[2] == 1:
                output.append("0N4R")
            # exon 2 present, exon 3 skipped, exon 10 present
            elif row[0] == 1 and row[1] == 0 and row[2] == 1:
                output.append("1N4R")
            # exon 2 skipped, exon 3 present, exon 10 present
            elif row[0] == 0 and row[1] == 1 and row[2] == 1:
                output.append("1N*4R")
            # exon 2 present, exon 3 present, exon 10 present
            elif row[0] == 1 and row[1] == 1 and row[2] == 1:
                output.append("2N4R")
            # exon 2 skipped, exon 3 skipped, exon 10 skipped
            elif row[0] == 0 and row[1] == 0 and row[2] == 0:
                output.append("0N3R")
            # exon 2 present, exon 3 skipped, exon 10 skipped
            elif row[0] == 1 and row[1] == 0 and row[2] == 0:
                output.append("1N3R")
            # exon 2 skipped, exon 3 present, exon 10 skipped
            elif row[0] == 0 and row[1] == 1 and row[2] == 0:
                output.append("1N*3R")
            # exon 2 present, exon 3 present, exon 10 skipped
            elif row[0] == 1 and row[1] == 1 and row[2] == 0:
                output.append("2N3R") 
            # exon 2 present, exon 3 skipped, no exon 10: truncated
            elif row[0] == 1 and row[1] == 0 and row[2] == 1001:
                output.append("1N_E10'") 
            elif row[0] == 0 and row[1] == 1 and row[2] == 1001:
                output.append("1N*_E10'") 
            # exon 2 skipped, exon 3 skipped, no exon 10: truncated
            elif row[0] == 0 and row[1] == 0 and row[2] == 1001:
                output.append("0N_E10'") 
            else:
                output.append(0)
    
    # Append the output to a new column and select relevant columns to return 
    mapt_exon_tab["MAPT_classification"] = output
    mapt_exon_tab = mapt_exon_tab[["MAPT_classification"]]
    
    # Tabulate the number of transcripts that fall into MAPT classifications
    mapt_exon_tab_counts = mapt_exon_tab.groupby(['MAPT_classification']).size()
    
    if "_" in '\t'.join(list(mapt_exon_tab.index)):
        print("Processing merged datasets")
        mapt_exon_tab["isoform"] = mapt_exon_tab.index
        mapt_exon_tab = reidentify_isoform_dataset(mapt_exon_tab)
        del mapt_exon_tab["isoform"]
        
    return(mapt_exon_tab, mapt_exon_tab_counts)
