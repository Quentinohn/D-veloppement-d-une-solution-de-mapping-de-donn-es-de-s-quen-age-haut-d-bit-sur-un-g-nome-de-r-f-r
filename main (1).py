from Bio import SeqIO
import cProfile
from tqdm import tqdm
import numpy as np
import json


def BWAo(str, rank, read,k,dico):
    index=len(read)-1
    s = read[index]
    band_start=dico[s][0][0]
    band_end=dico[s][0][-1]+1
    numDic = {"A":0, "C":1, "G":2, "T":3}
    for i in range(index-1,-1,-1):
        if(band_start == band_end):
            pass
        s = read[i]

        letterNum = numDic[s]
        if band_start>0:
            rank_top, rank_bottom = rank[band_start-1][letterNum], rank[band_end-1][letterNum]
        else:
            rank_top, rank_bottom = rank[band_start][letterNum], rank[band_end-1][letterNum]
        a=dico[s][0][0]
        band_start=a+rank_top
        band_end=a+rank_bottom
    return band_start,band_end


def kmers(str,k):
    l=len(str)
    b=l-l%k
    a=[]
    for i in range(0,b,k):
        a.append(str[i:i+k])
    if b!=l:
        a.append(str[b:])
    return a

def compinv(str):
    inv=str[::-1]
    mytable = inv.maketrans("ATGC", "TACG")
    return inv.translate(mytable)

def rank(str,char,i,j):
    Count1=0
    k=0
    for k in range(i):
        if str[k]==char:
            Count1+=1
    Count2 = Count1
    for k in range(i,j):
        if str[k]==char:
            Count2+=1

    return Count2, Count1


def band(str,char):
    a=np.array(sorted(str.strip()))
    return np.where(a==char)


def mapping(pos,S,k,kread):
    res=[]
    for j in range(len(pos)):
        err=0
        #errorlist=[]
        d=0
        for c in range(1,len(kread)):
            a,b=pos[j]+c*k+d,pos[j]+(c+1)*k+d

            if S[a:b]!=kread[c]:
                err-=1
                if err==-1:
                    break
                if c<len(kread)-1:
                    if S[a+k+1:b+k+1]==kread[c+1]:
                        d+=1
                    if S[a+k-1:b+k-1]==kread[c+1]:
                        d-=1
                #errorlist.append(a)
        if err!=-1:
            res.append(pos[j])
    return res

from Bio import SeqIO
it = SeqIO.parse('/home/azerty/Documents/boa/single.fq', 'fastq')
listesequenceread=[]
while True:
    try:
        seqRecord = next(it)
        listesequenceread.append(str(seqRecord.seq))
    except StopIteration:
        break
test=listesequenceread[0:10]

genome=[]
i=0
for seq_record in SeqIO.parse("/home/azerty/Documents/boa/genome.fna", "fasta"):
    a={}
    i+=1
    a[i]=str(seq_record.seq)
    genome.append(str(seq_record.seq))



file_reads = open("read.json")
reads=json.load(file_reads)

file_suffix_array = open("chr0.json",)
suffix_array=json.load(file_suffix_array)

file_bwt = open("Bwt8fin.json",)
bwt=json.load(file_bwt)

file_chrom = open("chromosome8.json")
gen=json.load(file_chrom)


Suffix_table,bwt, Ranks = d, file_bwt[0], file_bwt[1]
S=gen[0]
bgA=band(S,"A")
bgT=band(S,"T")
bgG=band(S,"G")
bgC=band(S,"C")

dictionnaire={"A":bgA,"T":bgT,"G":bgG,"C":bgC}


mapping=[]
for i in range(800000,900000):
    read=reads[i]
    k=5
    kread = kmers(read, k)
    align=BWAo(bwt, Ranks, kread[0],k,dictionnaire)
    pos=Suffix_table[align[0]:align[1]]
    map_result=mapping(pos,x[0],k,kread)
    if len(map_result)>0:
        mapping.append(map_result)

fichier=open("chromosome8mapping.json","w")
json.dump(mapping,fichier)
