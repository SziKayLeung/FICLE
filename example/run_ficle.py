#!/usr/bin/env python
# Szi Kay Leung (sl693@exeter.ac.uk)

"""
params:
    input_bed = bed file of final set of isoforms sorted and coloured 
    ref_dir = reference directory containing reference gtf files pertaining to each target gene of interest 
    input_gtf = gtf file of final set of isoforms (from SQANTI)
    noISM_path = SQANTI classification file with 3'ISM isoforms removed (as considered artefacts)
    ORF_dir = directory containing <prefix>_cpat.ORF_prob.best.tsv 
"""

##------------- Load modules -----------

import gtfparse
from gtfparse import read_gtf
import pandas as pd
import numpy as np
import collections
from collections import Counter
import csv
import sys
import dill
import sys
import shutil
import os
import argparse

srcPath = os.path.dirname(os.path.realpath(__file__)) +"/../src/"
sys.path.insert(1, srcPath)
import prepare_and_parse as prep
import matching_exons as me
import identify_exon_skipping as es
import alternative_promoter_terminator as apat
import identify_novel_exons as ne
import identify_intron_retention as ir
import identify_alternativeprime as aprime
import finalise_output as fo
import generate_multiregion as gm
import classify_mapt_isoforms as mapt
pd.options.mode.chained_assignment = None  # default='warn'


##------------- Input parameters -----------

gene = "RHBDF2"
root_dir = os.path.dirname(os.path.realpath(__file__)) 
i_bed = root_dir + "/input/input_sorted_bed.bed12"
i_gtf = root_dir + "/input/input_gtf.gtf"
noISM_path = root_dir + "/input/input_classification.txt"
gencode_gtf = root_dir + "/input/" + gene + "_gencode.gtf"
output_dir = root_dir + "/results/" + gene
output_log = root_dir + "/results/" + gene + "_parsed_transcripts.txt" 
output_log = open(output_log, "w")
orf_file = root_dir + "/input/input_cpat.ORF_prob.best.tsv"


##------------- Run FICLE -----------
print("************ Running FICLE...", file=sys.stdout)

# check the species dataset based on the gene input
# Human = Gene is all in capitals
if gene.isupper():
  species = "human"
else:
  species = "mouse"
print("Working with", species, "dataset")

# ORF for NMD Prediction
if orf_file:
  ORF = pd.read_csv(orf_file, sep= "\t")
else:
  print("Not using ORF for classification")
  ORF = None


# No ISM path 
if noISM_path:
  # keep only those that are retained from ISM filtering (i.e. removed 3'ISM from SQANTI)
  print("Keeping isoforms retained ISM filtering")
  noISM = pd.read_csv(noISM_path, sep = "\t")
else:
  noISM = None


# prepare directory
print("(Re)generating files")
fo.delete_make_dir(output_dir)


## read in gencode and transcriptome gtf
#gencode = pd.read_csv(ref+ gene + "_Manual_gencode.csv") 
gencode = prep.parse_gencode_reference(gencode_gtf, gene)
gencode.rename({'transcript': 'transcript_id'}, axis=1, inplace=True)

df = prep.parse_transcriptome_gtf(i_gtf, gene, noISM) 

# Parse through the transcriptome, classify and filter 
All_FilteredParsed = []
for count, transcript in enumerate(df['transcript_id'].unique()):
    #print(transcript)
    parsed = prep.parse_transcript(gencode,df,transcript,10, 100)
    All_FilteredParsed.append(prep.filter_parsed_transcript(gene,gencode,parsed, output_log))
    if count%50==0:
      print("Parsing through transcript", count)

# Aggregate all the filtered parsed output from each transcript into one big list
All_FilteredParsed = [x for l in All_FilteredParsed for x in l]

# QC: Check that the transcripts in the original transcriptome gtf captured in the final big list
if set(df["transcript_id"].unique()) != set(set([i.split(';',3)[0] for i in All_FilteredParsed])):
  print("Not all transcripts are annotated. Abort!.", file=sys.stderr)
  sys.exit(-1)

#All_Detected_Exons_Matching = identify_all_matching_exons(gencode, df, All_FilteredParsed)
AllKnownMatch, Mis2Match, MisMatch, MisjumpMatch, OtherClass, OtherClassMis2Match, SomeMatch, MisjumpMatch_NotAll = me.identify_all_matching_exons(gencode, df, All_FilteredParsed)
#All_Detected_Exons = identify_all_exons(gencode, df, All_FilteredParsed)
#All_Detected_Exons_NonMatching = list(set(All_Detected_Exons) - set(All_Detected_Exons_Matching))

# Tabulate exon presence 
print("Tabulating exon presence")
exon_tab = es.tabulate_exon_presence(gencode, df, All_FilteredParsed)

if gene == "MAPT" or gene == "Mapt":
  print("Further classifiying MAPT isoforms by exons 2, 3 and 10")
  mapt_exon_tab, mapt_exon_tab_counts = mapt.classify_mapt_isoforms(species, exon_tab, gencode)
  mapt_exon_tab.to_csv(output_dir + "/Stats/" + gene + "_further_classifications.csv",index_label="isoform")
  mapt_exon_tab_counts.to_csv(output_dir + "/Stats/" + gene + "_further_classifications_counts.csv")

# Exon Skipping 
print("Processing transcripts for exon skipping")
ES = es.identify_exon_skipping(gencode,exon_tab)
#if gene == "APOE": ES = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)
#if gene == "ABCA7": ES = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)
#if gene == "CD33": ES = gene_specific(gene, df, ES, All_FilteredParsed, gencode_gtf)
#ES = es.gene_specific_human(gene, df, ES, All_FilteredParsed, gencode_gtf)
ES = es.skip_not_AFexons(ES, gencode)

ES_Count, ES_SpecificExonSkipped, ES_Transcripts = es.output_exon_skipping_stats(ES)

# Alternative First Promoter    
Alternative_First_Promoter = apat.identify_alternative_promoter(df, ES, gene, gencode, gencode_gtf, All_FilteredParsed)


# Alternative Termination 
Alternative_Termination = apat.identify_alternative_termination(df, gencode, All_FilteredParsed)

# Alternative First Exon 
Alternative_First_exon = aprime.identify_alternative_first(df, All_FilteredParsed)

# Novel Exons 
print("Identifying transcripts with novel exons")
NE, NE_novel_co = ne.identify_novel_exon(df, gencode, All_FilteredParsed)
NE_classify, NExons_BeyondFirst, NExons_BeyondFirstLast, NExons_BeyondLast, NExons_Internal = ne.classify_novel_exon(gencode, df, NE)
NE_pertrans_classify_counts = pd.DataFrame()
if len(NE) > 0: 
  NE_pertrans_counts, NE_classify_counts, NE_pertrans_classify_counts = ne.novel_exon_stats(NE, NE_classify)
NE_classify_counts.to_csv(output_dir + "/Stats/" + gene + "_NE_counts.csv")
NE_pertrans_classify_counts.to_csv(output_dir + "/Stats/" + gene + "_NE_counts_pertrans.csv")
NE_novel_co.to_csv(output_dir + "/Stats/" + gene + "_NE_coordinates.csv")

# Intron Retention 
print("Identifying transcripts with intron retention")
IR, IR_Counts, IR_Transcripts, IR_Exon1, IR_LastExon = ir.identify_intron_retention(df, All_FilteredParsed, gencode)

# Alternative A5' and A3' 
print("Identifying transcripts with alternative 5' and 3' sites")
A5A3, A5A3_pertrans_counts, A5A3_Counts, A5A3_Transcripts = aprime.identify_A5A3_transcripts(df, All_FilteredParsed)

# Final Output 
# Event lists = each category is generated from previous functions and contain list of transcripts under that event type
categories = [AllKnownMatch, Mis2Match, MisMatch, MisjumpMatch, OtherClass, OtherClassMis2Match, SomeMatch,MisjumpMatch_NotAll,
              A5A3_Transcripts, Alternative_First_exon, Alternative_First_Promoter, Alternative_Termination,
              ES_Transcripts, IR_Transcripts, IR_Exon1, IR_LastExon, NExons_BeyondFirst, NExons_Internal,NExons_BeyondLast,NExons_BeyondFirstLast]

Transcript_Classifications = fo.generate_aggregated_classification(df, categories)
fo.prioritise_write_output(df, Transcript_Classifications, output_dir, gene, i_bed)
Transcript_Classifications_Remapped = fo.populate_classification(Transcript_Classifications, A5A3, IR_Counts, ES_Count, NE, NE_pertrans_classify_counts)

#NMD_output = call_NMD_prediction()

# generate multiregion file 
#gm.generate_multiregion(df, NE_novel_co,gencode)

# All other stats 
gencode.to_csv(output_dir + "_gencode_automated.csv")
exon_tab.to_csv(output_dir + "/Stats/" + gene + "_Exon_tab.csv")
ES.to_csv(output_dir + "/Stats/" + gene + "_Exonskipping_generaltab.csv")
A5A3.to_csv(output_dir + "/Stats/" + gene + "_A5A3_tab.csv")
IR.to_csv(output_dir + "/Stats/" + gene + "_IntronRetention_tab.csv", index = False)
ES_SpecificExonSkipped.to_csv(output_dir + "/Stats/" + gene + "_Exonskipping_tab.csv")
#NMD_output.to_csv(output_dir + "/Stats/" + gene + "_NMD_orf.csv")
Transcript_Classifications_Remapped.to_csv(output_dir + "/Stats/" + gene + "_Final_Transcript_Classifications.csv")
gencode.to_csv(output_dir + "/Stats/" + gene + "flattened_gencode.csv")

print("**** All Done! ****")
