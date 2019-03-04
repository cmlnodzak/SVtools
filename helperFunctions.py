#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 15:51:52 2019

@author: cmlnodzak
"""
from math import floor, log10, fabs
from pybedtools import BedTool
import pandas as pd
import re


# return the order of magnitude of a number.
def orderOfMag(val):
    n = floor(log10(fabs(val)))
    return n

# URL = ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz
    
#gtf ='gencode.v25.annotation.gtf'

# extract pandas dataframe from a GTF for a particular feature of interest.

def _attributeExtractor(attr_df):
    '''
    =================================================================+========
    Private method used by gtfReader() method to extract attributes 
    from a Pandas DataFrame. Returns a 12 column DataFrame with original 
    'attributes' column dropped.
    ==========================================================================
    '''
    attr_df['gene_id'] = attr_df['attributes'].map(lambda x: re.search(r'gene_id "ENSG\d+.\d+', x)[0].split()[1].strip("\""))    
    attr_df['gene_type'] = attr_df['attributes'].map(lambda x: re.search(r'gene_type\s"\w+',x)[0].split()[1].strip("\""))
    attr_df['gene_name'] = attr_df['attributes'].map(lambda x: re.search(r'gene_name\s"\w+',x)[0].split()[1].strip("\""))
    attr_df['gene_status'] = attr_df['attributes'].map(lambda x: re.search(r'gene_status\s"\w+',x)[0].split()[1].strip("\""))
    attr_df = attr_df.drop('attributes',axis=1)
    return attr_df

def gtfReader(gtf, feat, group = None, chrom = None):
    '''
    ==========================================================================
    gtfReader: a function to parse gtf files into clean pandas Dataframes.
    ==========================================================================
    parameters:
    ==========================================================================
        -- gtf: a gene transfer format file. 
            (see here ensembl.org/info/website/upload/gff.html for more.)
        -- feat: a feature of interest, may be one of 'gene', 'transcript', 
                'exon'.
        -- group: optional, the gene_type that the feature maps to. may be 
            'protein_coding, lncRA', etc. 
            For full list run gtfReader() then df.gene_type.drop_duplicates().
        -- chrom: optional, the name of the chromosome found in the '
                seqname' field.
    ==========================================================================
    output:
    ==========================================================================
        Returns a Pandas DataFrame with column names listed below.
        column names = ['seqname', 'start', 'end', 
                            'score', 'strand', 'frame', 'gene_id', 
                            'gene_type', 'gene_name', 'gene_status']
    ==========================================================================
    '''
    anno_df = BedTool(gtf).to_dataframe(skiprows=5)
    anno_df = anno_df[anno_df.feature == feat]
    anno_df = _attributeExtractor(anno_df)
    if group == 'lncRNA':
        group = ['lincRNA', 'antisense', 'sense_overlapping', 'macro_lncRNA', 'sense_intronic', 'bidirectional_promoter_lncRNA' ]
        anno_df = anno_df[anno_df['gene_type'].isin(group)].reset_index(drop=True)   
    elif group != None:
        anno_df = anno_df[anno_df.gene_type == group].reset_index(drop=True)
    if chrom != None:
        anno_df[anno_df.seqname == chrom].reset_index(drop=True)
    return anno_df.drop(['source', 'feature'],axis=1)


class VCFparser:
    col_names =  ['#chrom', 'start', 'stop', 'REF', 'ALT', 'QUAL', 'FILT', 'sample', 'GT']
    hgsv_list = ['HG00512','HG00513','HG00514', 'HG00731', 'HG00732', 'HG00733', 'NA19238', 'NA19239', 'NA19240']

    def __init__(self,infile, outfy, genoFilt=None):
        self.genoFilt=genoFilt
        self.bedtoolsList = []
        with open(infile) as SVs, open(outfy,'w') as outfile:
                data = []
                for line in SVs.readlines():
                        if line.startswith('#'):
                                continue
                        else:
                                line = line.rstrip().split('\t')
                                end = line[7].split(';')[1].split('=')[1]
                                line.pop(2)
                                line.insert(2,end)
                                info = line[7].split(';')[7].split(',')
                                for genos in info:
                                        geno = genos.split(':')
                                        sampid = str(geno[4])
                                        genotype = str(geno[3])
                                        SVdata=line[0:7]
                                        SVdata.append(sampid)
                                        SVdata.append(genotype)
                                        data.append(SVdata)
                                        keep="\t".join(SVdata)
                                        outfile.write(keep+'\n')                  
                self.frame = pd.DataFrame.from_records(data)
                self.frame.columns = self.col_names
        self.bedtoolsList = self._Frame2Bed(self.frame,self.hgsv_list,self.genoFilt)
   
     
#pass in a list of sample names as they are written in the file
    def _Frame2Bed(self,df,samples,genoFilt):
        bedList = []
        if self.genoFilt == "heterozygous":
        # extract heterozygous deletions and insertions, group by samples.
            df = df.loc[df['GT'] == '0/1']
        df = df.groupby('sample')
        for i in samples:
            i = df.get_group(i).drop_duplicates()
            i_bt = BedTool.from_dataframe(i)
            bedList.append(i_bt)
        return bedList

