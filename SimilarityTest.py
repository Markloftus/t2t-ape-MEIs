import pandas as pd
import os
import numpy as np
import ast
import pysam
import collections
from Bio.Seq import Seq

gorillaDF = pd.read_csv("/project/mkonkel/tangeno/TEleaves/LS_Primates/dataframes/hs1_Matches_06-12-2024.csv")
GenomeDict2={
'GorGor-mat':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028885495.2_NHGRI_mGorGor1-v2.0_mat_genomic.fna',
    'GorGor-pat':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028885475.2_NHGRI_mGorGor1-v2.0_pat_genomic.fna',
    'PanPan-mat':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028858845.2_NHGRI_mPanPan1-v2.0_mat_genomic.fna',
    'PanPan-pat':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028858825.2_NHGRI_mPanPan1-v2.0_pat_genomic.fna',
    'PanTro-pri':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna',
    'PanTro-alt':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028858805.2_NHGRI_mPanTro3-v2.0_alt_genomic.fna',
    'PonAbe-pri':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna',
    'PonAbe-alt':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028885685.2_NHGRI_mPonAbe1-v2.0_alt_genomic.fna',
    'PonPyg-alt':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028885525.2_NHGRI_mPonPyg2-v2.0_alt_genomic.fna',
    'PonPyg-pri':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna',
    'SymSyn-alt':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028878085.2_NHGRI_mSymSyn1-v2.0_alt_genomic.fna',
    'SymSyn-pri':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028878055.2_NHGRI_mSymSyn1-v2.0_pri_genomic.fna',
    'hs1.txt':'/home/loftus/genomes/chm13v2/chm13v2.fasta'
}
genomeDict ={
 'GorGor':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna',
 'PanPan':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna',
 'PanTro':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna',
 'PonAbe':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna',
 'PonPyg':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna',
 'SymSyn':'/project/mkonkel/tangeno/TEleaves/genomes/2024_Genomes/GCA_028878055.2_NHGRI_mSymSyn1-v2.0_pri_genomic.fna',
 'hs1':'/home/loftus/genomes/chm13v2/chm13v2.fasta',   
}
T2TDict={'NC_060925.1':'chr1',
'NC_060926.1':'chr2',
'NC_060927.1':'chr3',
'NC_060928.1':'chr4',
'NC_060929.1':'chr5',
'NC_060930.1':'chr6',
'NC_060931.1':'chr7',
'NC_060932.1':'chr8',
'NC_060933.1':'chr9',
'NC_060934.1':'chr10',
'NC_060935.1':'chr11',
'NC_060936.1':'chr12',
'NC_060937.1':'chr13',
'NC_060938.1':'chr14',
'NC_060939.1':'chr15',
'NC_060940.1':'chr16',
'NC_060941.1':'chr17',
'NC_060942.1':'chr18',
'NC_060943.1':'chr19',
'NC_060944.1':'chr20',
'NC_060945.1':'chr21',
'NC_060946.1':'chr22',
'NC_060947.1':'chrX',
'NC_060948.1':'chrY'}
T2TDict2 = {y:x for x,y in T2TDict.items()}
seqSims=[]
for row in gorillaDF.index:
    mykmerSize=14
    similarities=[]
    
    contig = T2TDict2[str(gorillaDF.at[row,'Loci'].split(":")[0])]
    
    start = str(int(gorillaDF.at[row,'Loci'].split(":")[1].split("-")[0])-500)
    end = str(int(gorillaDF.at[row,'Loci'].split(":")[1].split("-")[1])+500)
    
    sequence = ''.join(pysam.faidx(genomeDict['hs1'], str(contig)+":"+str(start)+"-"+str(end)).split()[1:])
    
    for match in ast.literal_eval(str(gorillaDF.at[row,'Matches'])):
        
        querySpecies = str(gorillaDF.at[row,'Query_Species'])
        
        if querySpecies == 'hs1.txt':
            coordinate = str('_'.join(match.split("_")[1:]))
        else:
            coordinate = str(match.split("_")[1])
        
        if str(match.split("_")[0]) == '+':
            
            sequence2 = ''.join(pysam.faidx(GenomeDict2[querySpecies], coordinate).split()[1:])
            
        else:
            presequence = Seq(''.join(pysam.faidx(GenomeDict2[querySpecies], coordinate).split()[1:]))
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
gorillaDF.to_csv("/project/mkonkel/tangeno/TEleaves/LS_Primates/dataframes/hs1_Matches_Filled_06-12-2024.csv")
