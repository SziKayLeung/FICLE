#!/usr/bin/env python3
# FICLE
# Authors: Szi Kay Leung

__author__  = "S.K.Leung@exeter.ac.uk"
__version__ = '1.1.3'  # Python 3.7

## Load Libraries
import gtfparse
from gtfparse import read_gtf
import pandas as pd
import numpy as np
import collections
from collections import Counter
import csv
import sys
import sys
import shutil
import os
import argparse

sys.path.append("/src")
from src import prepare_and_parse as prep
from src import matching_exons as me
from src import identify_exon_skipping as es
from src import alternative_promoter_terminator as apat
from src import identify_novel_exons as ne
from src import identify_intron_retention as ir
from src import identify_alternativeprime as aprime
from src import finalise_output as fo
from src import generate_multiregion as gm
from src import classify_mapt_isoforms as mapt
pd.options.mode.chained_assignment = None  # default='warn'

def annotate_gene(args):
    
    # prepare output directories
    args = prep.create_output_dir(args)
    
    # subset reference from gene
    prep.subset_reference(args)
    output_log = open(args.gene_stats_dir + args.genename + "parsed_transcripts.txt", "w") 
    
    # check the species dataset based on the gene input
    # Human = Gene is all in capitals
    if args.genename.isupper():
        species = "human"
    else:
        species = "mouse"
    print("Working with", species, "dataset")

    # ORF for NMD Prediction
    if args.cpat:
        ORF = pd.read_csv(args.cpat, sep= "\t")
    else:
        print("Not using ORF for classification")
        ORF = None

    # read in gencode and transcriptome gtf
    gencode,order = prep.parse_gencode_reference(args)
    gencode.rename({'transcript': 'transcript_id'}, axis=1, inplace=True)

    df = prep.parse_transcriptome_gtf(args, order) 
    
    if len(df) > 0:

        # Parse through the transcriptome, classify and filter 
        All_FilteredParsed = []
        for count, transcript in enumerate(df['transcript_id'].unique()):
            #print(transcript)
            parsed = prep.parse_transcript(gencode,df,transcript,10, order)
            All_FilteredParsed.append(prep.filter_parsed_transcript(gencode,parsed, output_log))
            if count%50==0:
                print("Parsing through transcript", count)

        # Aggregate all the filtered parsed output from each transcript into one big list
        All_FilteredParsed = [x for l in All_FilteredParsed for x in l]

        # QC: Check that the transcripts in the original transcriptome gtf captured in the final big list
#        if set(df["transcript_id"].unique()) != set(set([i.split(';',3)[0] for i in All_FilteredParsed])):
        if np.all(np.isin([i.split(';',3)[0] for i in All_FilteredParsed], df["transcript_id"].unique()))==False:
            print("Mismatch transcripts")
            print(set(df["transcript_id"].unique()))
            print(set(set([i.split(';',3)[0] for i in All_FilteredParsed])))
            sys.exit(-1)

	df=df[df["transcript_id"].isin([i.split(';',3)[0] for i in All_FilteredParsed])]

        # Check for matching or somematching 
        AllKnownMatch, SomeMatch = me.identify_all_matching_exons(gencode, df, All_FilteredParsed)

        # Tabulate exon presence 
        print("Tabulating exon presence")
        exon_tab = es.tabulate_exon_presence(args, gencode, df, All_FilteredParsed)

        if args.genename == "MAPT" or args.genename == "Mapt":
            print("Further classifiying MAPT isoforms by exons 2, 3 and 10")
            mapt_exon_tab, mapt_exon_tab_counts = mapt.classify_mapt_isoforms(species, exon_tab, gencode)
            mapt_exon_tab.to_csv(args.gene_stats_dir + args.genename + "_further_classifications.csv",index_label="isoform")
            mapt_exon_tab_counts.to_csv(args.gene_stats_dir + args.genename + "_further_classifications_counts.csv")

        # Exon Skipping 
        print("Processing transcripts for exon skipping")
        ES = es.identify_exon_skipping(gencode,exon_tab,All_FilteredParsed)
        ES = es.skip_not_AFexons(ES, gencode)
        ES_Counts, ES_Transcripts = es.output_exon_skipping_stats(args, ES)

        # Alternative First Promoter    
        Alternative_First_Promoter = apat.identify_alternative_promoter(df, gencode, All_FilteredParsed)

        # Alternative Termination 
        Alternative_Termination = apat.identify_alternative_termination(df, gencode, All_FilteredParsed)

        # Alternative First Exon 
        AF, AF_Counts, AF_Transcripts = aprime.identify_alternative_first(df, All_FilteredParsed)

        # Novel Exons 
        print("Identifying transcripts with novel exons")
        NE = ne.identify_novel_exon(args, df, gencode, All_FilteredParsed)
        NE_classify, NExons_BeyondFirst, NExons_BeyondFirstLast, NExons_BeyondLast, NExons_Internal, NExons_First, NExons_Last = ne.classify_novel_exon(gencode, order, df, NE, All_FilteredParsed)
        NE_Counts = pd.DataFrame()
        if len(NE) > 0: NE_Counts = ne.novel_exon_stats(args, NE, NE_classify)

        # Intron Retention 
        print("Identifying transcripts with intron retention")
        IR_Counts, IR_Transcripts, IR_Exon1, IR_LastExon = ir.identify_intron_retention(args, df, All_FilteredParsed, gencode)

        # Alternative A5' and A3' 
        print("Identifying transcripts with alternative 5' and 3' sites")
        A5A3_Counts, A5A3_Transcripts = aprime.identify_A5A3_transcripts(args, df, All_FilteredParsed)

        # Final Output 
        # Event lists = each category is generated from previous functions and contain list of transcripts under that event type
        categories = [AllKnownMatch, SomeMatch, A5A3_Transcripts, AF_Transcripts, 
              Alternative_First_Promoter, Alternative_Termination,
              ES_Transcripts, IR_Transcripts, IR_Exon1, IR_LastExon, NExons_BeyondFirst,
              NExons_Internal,NExons_BeyondLast,NExons_BeyondFirstLast]

        Transcript_Classifications = fo.generate_aggregated_classification(df, categories)
        fo.prioritise_write_output(args, df, Transcript_Classifications)
        fo.populate_classification(args, Transcript_Classifications, A5A3_Counts, IR_Counts, ES_Counts, NE_Counts)   
    
    else:
        print("No detected isoforms")
        shutil.rmtree(args.gene_root_dir)
    
    print("**** All Done! ****")

    
def main():
    parser = argparse.ArgumentParser(description="Full Isoform Characterisation from (Targeted) Long-read Experiments")
    parser.add_argument("-n","--genename", help='\t\tTarget gene symbol')
    parser.add_argument("-r","--reference", help='\t\tGene reference annotation (<gene>_gencode.gtf)')
    parser.add_argument("-b","--input_bed", help='\t\tInput bed file of all the final transcripts in long-read derived transcriptome.')
    parser.add_argument("-g","--input_gtf", help='\t\tInput gtf file of all the final transcripts in long-read derived transcriptome.')
    parser.add_argument("-c","--input_class", help='\t\tSQANTI classification file')
    parser.add_argument("--cpat", help='\t\ORF_prob.best.tsv file generated from CPAT',required=False)
    parser.add_argument("-o","--output_dir", help='\t\tOutput path for the annotation and associated files')
    parser.add_argument("-v","--version", help="Display program version number.", action='version', version='FICLE '+str(__version__))
    
    args = parser.parse_args()
    print("************ Running FICLE...", file=sys.stdout)
    print("version:", __version__)
    annotate_gene(args)        
    
if __name__ == "__main__":
    main()
