def extract_co(df, gencode, num):
    chrom = df.loc[df["transcript_id"] == transcript,"seqname"].values[0]
    dat = gencode.loc[gencode["updated_exon_number"] == num,]
    dat["diff"] = abs(dat["start"] - dat["end"])
    start = dat[dat['diff']==dat['diff'].max()]["start"].values[0]
    end = dat[dat['diff']==dat['diff'].max()]["end"].values[0]
    return([chrom, start, end, num])

def extract_co_novel(df, NE_novel):
    chrom = df.loc[df["transcript_id"] == transcript,"seqname"].values[0]
    start = int(NE_novel.split(",",2)[0])
    end = int(NE_novel.split(",",2)[1]) 
    diff = abs(start - end)
    return([chrom, start, end, diff])

def generate_multiregion(df, NE_novel_co, gencode):
    novel_output = []
    for i in NE_novel_co:
        novel_output.append(extract_co_novel(df,gencode, i))
    
    output = []
    for i in range(1, max([int(i) for i in gencode["updated_exon_number"]])):
        output.append(extract_co(df, gencode, str(i)))
    
    final = pd.concat([pd.DataFrame(output), pd.DataFrame(novel_output)], axis=0)
    final = final.sort_values(by=[1])
    final.to_csv(output_dir + "_multiregon.bed", sep = "\t", index = False, header = False)