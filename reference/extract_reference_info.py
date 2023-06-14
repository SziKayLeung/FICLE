#!/usr/bin/env python3
# Szi Kay Leung (sl693@exeter.ac.uk)

"""
Aim: Split reference genome annotation gtf by target gene name 
"""

## Load Libraries
import argparse
from gtfparse import read_gtf
import pandas as pd
import csv
import numpy as np

# switch off loc warnings
pd.options.mode.chained_assignment = None  # default='warn'


def writeGTF(inGTF,file_path):
    """
    Write a GTF dataframe into a file
    :param inGTF: GTF dataframe to be written. It should either have 9 columns with the last one being the "attributes" section or more than 9 columns where all columns after the 8th will be colapsed into one.
    :param file_path: path/to/the/file.gtf
    :returns: nothing
    https://github.com/mpg-age-bioinformatics/AGEpy/blob/master/AGEpy/gtf.py
    """
    cols=inGTF.columns.tolist()
    if len(cols) == 9:
        if 'attribute' in cols:
            df=inGTF
    else:
        df=inGTF[cols[:8]]
        df['attribute']=""
        for c in cols[8:]:
            if c == cols[len(cols)-1]:
                df['attribute']=df['attribute']+c+' "'+inGTF[c].astype(str)+'";'
            else:
                df['attribute']=df['attribute']+c+' "'+inGTF[c].astype(str)+'"; '
    df.to_csv(file_path, sep="\t",header=None,index=None,quoting=csv.QUOTE_NONE)

    
def subsetGTF(gencode_gtf, gene, output_dir):
    print("Filtering for:", gene)
    filtered_gtf = gencode_gtf.loc[gencode_gtf["gene_name"] == gene,]
    output_path = output_dir + "/" + gene + "_gencode.gtf"
    writeGTF(filtered_gtf,output_path)    
    

def find_gene_length_exon_num(gencode_gtf, gene):
    '''
    The gene length, longest transcript length, minimum and maximum   
    '''

    print("Extracting information about", gene)
    # subset gencode by gene
    gencode = gencode_gtf.loc[gencode_gtf['gene_name'] == gene,]
    gencode.replace("", str("NaN"), inplace=True)
    
    # Description Feature determination 
    # length of all features including exons, transcript 
    # iso = number of isoforms based on the transcript id
    # number of exons for that gene
    gencode["length"] = abs(gencode["start"] - gencode["end"])
    iso = [i for i in list(gencode["transcript_id"].unique()) if i != "NaN"]
    exon = [int(i) for i in list(gencode["exon_number"].unique()) if i != "NaN"]
    gene_length = gencode[gencode["feature"] == "transcript"]["length"]
    
    # length of each transcript calculated from sum of the exon length across each transcript
    length = []
    for i in iso:
        length.append(sum(gencode[(gencode["feature"] == "exon") & (gencode["transcript_id"] == i)]["length"]))

    # function for output of minimum and maximum value 
    def output(val):
        return(str(min(val)) + "-" + str(max(val)))        

    return ([gene, len(iso), output(length), output(exon), max(gene_length), max(length), max(exon)])




def identify_constitutive_alternative_exon(args, gene):
    
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
    #rnaseq_gene_normalised_counts, excluded_transcripts, monoexons, gencode_dir = input_files(dataset)
    
    if args.exc is not None:
        print("**** Removing mono-exonic transcripts.")
        excluded_transcript = pd.read_csv(args.exc)["excluded_transcript"].values
    else:
        excluded_transcript = []
        
    
    # Read in the gencode classification (previously manually curated to flatten the gene structure)
    gencode = pd.read_csv(args.split + "/" + gene + "_gencode_automated.csv") 
    gencode["start_end"] = gencode["start"].astype(str) + "_" + gencode["end"].astype(str)
    #print(gencode["transcript"].unique())
    
    # Remove monoexons to not loop through the unique exons for complication
    #monoexons = ["ENSMUST00000207525.1","ENSMUST00000207106.1","ENSMUST00000207596.1",'ENSMUST00000190478.1'] # Apoe, Picalm
    #gencode = gencode.loc[~gencode["transcript"].isin(monoexons),]
    
    # Loop through each set of start and end coordinates
    # Then loop through each set of transcripts start and end coordinates
    alternative_exon = []
    for e in gencode["start_end"]:
        e_start = e.split("_")[0]
        e_end = e.split("_")[1]
        
        # loop through each transcript 
        for t in gencode["transcript_id"].unique(): 
            if t not in excluded_transcript: # not to include the short transcripts 
                trans_exon = gencode.loc[gencode["transcript_id"] == t,"start_end"].values
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


def gene_description(args, gencode_gtf):
    '''
    run each gene through to:
    1/ find_gene_length_exon_num and merge with the counts results for expression
    2/ identify the constitutve and alternative exon
    '''

    # loop through each target gene to dataframe
    output_length = []
    output_exon = [] 
    for g in args.glist:
        output_length.append(find_gene_length_exon_num(gencode_gtf, g))
        output_exon.extend(identify_constitutive_alternative_exon(args, g))
    ### Length and expression
    target_des = pd.DataFrame(output_length, columns = ["associated_gene","Number_of_Gencode_Isoforms","Transcript_Length","Number_of_Exons","MaxGeneLength","MaxTransLength","Maxexons"])
    print(target_des)
    
    # read the normalised counts and determine average for each target gene
    #rnaseq_gene_normalised_counts, excluded_transcripts, monoexons, gencode_dir = input_files(dataset)
    if args.short_read is not None:
        print("**** Reading short-read normalised counts.")
        normalised_counts = pd.read_csv(args.short_read,index_col=0)
        normalised_counts['MeanRNASeqCounts'] = normalised_counts.mean(axis=1)
        normalised_counts["MeanRNASeqCounts"] = round(normalised_counts["MeanRNASeqCounts"],2) # round counts to 2 dp
        normalised_counts = normalised_counts.loc[normalised_counts.index.isin(args.glist)]["MeanRNASeqCounts"]
    else:
        normalised_counts = pd.DataFrame()
    
    # merge results of length (keep on the target_des even if no normalised counts)
    final_length = pd.merge(target_des, normalised_counts, left_on = "associated_gene", right_index=True, how='left')
    
    ### Constitutive and alternative exons 
    final_exon = pd.DataFrame(output_exon, columns = ["exon","exon_status","associated_gene"])
    return final_length, final_exon


def main():
    parser = argparse.ArgumentParser(description="Splitting gencode reference by gene for easier handle")
    parser.add_argument('--r', "--reference_gtf", help='\t\tPath to Gencode Reference gtf.')
    parser.add_argument("--glist", nargs="+", required=True, help='\t\tList of Target Genes')
    parser.add_argument('--g', "--gene", required=False, help='\t\tTarget gene for filtering gencode reference.')
    #parser.add_argument('--d', "--dataset", help='\t\tDataset essential for using the right input files: Mouse or ADBDR.')
    parser.add_argument("--split", help='\t\tPath of directory containing manual.csv')
    parser.add_argument("--exc", required=False,help='\t\tPath of fiile containing mono-exonic isoforms to exclude')
    parser.add_argument('--short_read', required=False, help='\t\tShort-read RNA-Seq normalised counts.')
    parser.add_argument('--o',"--output_dir", default=None, required=False, help='\t\tOutput directory, default: split directory')
    
    args = parser.parse_args()
    print("Reading in:", args.r)
    gencode_gtf = read_gtf(args.r)
    
    if args.o is None:
        args.o = args.split
    
    print("Filtering reference genome")
    if args.glist:
        Final_Stats = []
        #for gene in args.glist: 
          #subsetGTF(gencode_gtf, gene, args.o)
      
        # Concatenate results       
        Final_length_df, Final_exon_df = gene_description(args, gencode_gtf)
        Final_length_df.to_csv(args.o + "/" + "TargetGene_Reference_LengthNum.csv")
        Final_exon_df.to_csv(args.o + "/" + "TargetGene_Reference_AltConExons.csv")
    else:
        subsetGTF(gencode_gtf, args.g, args.o)

if __name__ == "__main__":
    main()
