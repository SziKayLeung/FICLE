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
    

def main():
    parser = argparse.ArgumentParser(description="Splitting gencode reference by gene for easier handle")
    parser.add_argument('--r', "--reference_gtf", help='\t\tPath to Gencode Reference gtf.')
    parser.add_argument("--glist", nargs="+", required=False, help='\t\tList of Target Genes')
    parser.add_argument('--g', "--gene", help='\t\tTarget gene for filtering gencode reference.')
    parser.add_argument('--o',"--output_dir", help='\t\tOutput directory')
   
    args = parser.parse_args()
    print("Reading in:", args.r)
    gencode_gtf = read_gtf(args.r)
    
    print("Filtering reference genome")
    if args.glist:
      Final_Stats = []
      for gene in args.glist: 
        subsetGTF(gencode_gtf, gene, args.o)
        
    else:
      subsetGTF(gencode_gtf, args.g, args.o)
    
if __name__ == "__main__":
    main()
