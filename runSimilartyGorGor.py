import pandas as pd
import os
import numpy as np
import ast
import pysam
import collections
from Bio.Seq import Seq

gorillaDF = pd.read_csv("/project/mkonkel/tangeno/TEleaves/LS_Primates/dataframes/GorGor_Matches.csv")

genomeDict ={
 'GorGor':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna',
 'PanPan':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna',
 'PanTro':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna',
 'PonAbe':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna',
 'PonPyg':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna',
 'SymSyn':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028878055.2_NHGRI_mSymSyn1-v2.0_pri_genomic.fna',
 'hs1':'/home/loftus/genomes/chm13v2/chm13v2.fasta',   
}

seqSims=[]
for row in gorillaDF.index:
    mykmerSize=14
    similarities=[]
    
    contig = str(gorillaDF.at[row,'Loci'].split(":")[0])
    
    start = str(int(gorillaDF.at[row,'Loci'].split(":")[1].split("-")[0])-500)
    end = str(int(gorillaDF.at[row,'Loci'].split(":")[1].split("-")[1])+500)
    
    sequence = ''.join(pysam.faidx(genomeDict['GorGor'], str(contig)+":"+str(start)+"-"+str(end)).split()[1:])
    
    for match in ast.literal_eval(str(gorillaDF.at[row,'Matches'])):
        
        querySpecies = str(gorillaDF.at[row,'Query_Species'])
        
        if querySpecies == 'hs1':
            coordinate = str('_'.join(match.split("_")[1:]))
        else:
            coordinate = str(match.split("_")[1])
        
        if str(match.split("_")[0]) == '+':
            
            sequence2 = ''.join(pysam.faidx(genomeDict[querySpecies], coordinate).split()[1:])
            
        else:
            presequence = Seq(''.join(pysam.faidx(genomeDict[querySpecies], coordinate).split()[1:]))
            sequence2 = str(presequence.reverse_complement())

            
        kmers=[]
        insertionkmers=[]
        i=0
        while i < (len(sequence)-(mykmerSize-1)):
            kmers.append(str(sequence[i:i+mykmerSize]).upper())
            insertionkmers.append(str(sequence2[i:i+mykmerSize]).upper())
            i+=1

        uniqueKmers = set(kmers+insertionkmers)
        tempDF = pd.DataFrame(0, index=['Reference','MEI'], columns=[x for x in uniqueKmers])

        for kmer in kmers:
            tempDF.at['Reference',kmer]+=1
        for kmer2 in insertionkmers:
            tempDF.at['MEI',kmer2]+=1

        tempDF.loc['Diff'] = tempDF.loc['Reference'] - tempDF.loc['MEI']
        tempDF.loc['Sum'] = tempDF.loc['Reference'] + tempDF.loc['MEI']
        similarities.append(sum(abs(tempDF.loc['Diff']))/sum(abs(tempDF.loc['Sum'])))
        
    seqSims.append(similarities)
gorillaDF['Hit_Similarity']=seqSims
gorillaDF.to_csv("/project/mkonkel/tangeno/TEleaves/LS_Primates/dataframes/GorGor_Matches_Filled.csv")