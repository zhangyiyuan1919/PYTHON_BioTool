#DNA工具包
import collections
from structures import *
from collections import Counter
import random

def validateSeq(dna_seq):# 工具1：定义合法DNA碱基（支持IUPAC简并碱基可自行扩展）
    """
    严格校验DNA序列：
    - 若包含任何非A/T/C/G的字符，直接抛出ValueError，明确提示非法字符
    - 全合法则返回大写的标准序列
    """
    tmpseq = dna_seq.upper() # 收集所有非法字符
    invalid_chars = {nuc
                     for nuc in tmpseq
                     if nuc not in Nucleotides
                     }
    if invalid_chars:
        raise ValueError(f"❌ DNA序列包含非法字符: {sorted(invalid_chars)}")
    return tmpseq

def countFrequency(seq):#工具2：统计碱基数量
    tmpFreqDict={"A":0,"T":0,"C":0,"G":0} #计数字典
    for nuc in seq:
        tmpFreqDict[nuc]+=1
    return tmpFreqDict
#另一种解法
#def countNucFrequency(seq):
#    return dict(collections.Counter(seq))
def transcription(seq):#工具3：序列转录 DNA→RNA
    """DNA→RNA转录，用U代替T"""  #在自定工具包中相当于加入注释
    return seq.replace("T","U")

def reverse_complement(seq):#工具4：互补序列
    """DNA的互补链，顺序是5‘→3’"""
   # return ''.join([DNA_ReverseComplement[Nuc] for Nuc in seq])[::-1]#[::-1]首位调转互补后的序列
    # PY方法更快捷一些
    mapping=str.maketrans('ATCG','TAGC')
    return seq.translate(mapping)[::-1]


#工具5：基因密度统计➡gc含量统计
def gc_content(seq):
    """在一段DNA或者RNA上GC含量"""
    return round((seq.count('C')+seq.count('G'))/len(seq)*100)
#滑动窗口 GC 含量统计
def gc_content_subsec(seq,k=20):
    res=[]
    for i in range(0,len(seq)-k+1,k):
        subseq=seq[i:i+k]
        res.append(gc_content(subseq))
    return res

#工具6：DNA密码子
def translate_seq(seq,init_pos=0):
    return [DNA_Codons[seq[pos:pos+3]] for pos in range(init_pos,len(seq)-2,3)]

#工具7：给出 DNA 序列中编码特定氨基酸的每个密码子的出现频率
def codon_usage(seq,aminoacid):

    tmpList=[] #创建空列表，用于存储所有符合条件的密码子
    for i in range(0,len(seq)-2,3): #步长为 3：按密码子（3 个碱基为 1 个单位）遍历 DNA 序列 确保最后一次取子串 seq[i:i+3] 不会越界
        if DNA_Codons[seq[i:i+3]]==aminoacid:
            tmpList.append(seq[i:i+3])

    freqDict=dict(Counter(tmpList))
    totalWight=sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq]=round(freqDict[seq]/totalWight,2)
    return  freqDict

#工具8：生成 DNA 序列的氨基酸序列，包括反向互补序列
def gen_reading_frames(seq):
    frames=[]
    frames.append(translate_seq(seq,0))
    frames.append(translate_seq(seq,1))
    frames.append(translate_seq(seq,2))
    #互补链的氨基酸序列
    frames.append(translate_seq(reverse_complement(seq), 0))
    frames.append(translate_seq(reverse_complement(seq), 1))
    frames.append(translate_seq(reverse_complement(seq), 2))
    return frames

#工具9：提取蛋白质序列 所有可能的序列
def proteins_from_rf(aa_seq): #从gen_reading_frames中找到所有基因序列
    current_prot=[]
    proteins=[]
    for aa in aa_seq:
        if aa =="_": #发现终止子代表的aa
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot=[]
        else: #发现起始密码子代表的M
            if aa=="M":
                current_prot.append("")
            for i in  range(len(current_prot)):
                current_prot[i]+=aa
    return proteins

#工具10：提取所有蛋白质序列
def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0, ordered=False):
    """计算所有开放阅读框对应的全部可能蛋白质"""
    """Protine Search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
    """API can be used to pull protein info"""
    if endReadPos > startReadPos:
        rfs=gen_reading_frames(seq[startRead:endRead])
    else:
        rfs=gen_reading_frames(seq)

    res=[]
    for rf in rfs:
        prots=proteins_from_rf(rf)
        for p in prots:
            res.append(p)

    if ordered:
        return sorted(res,key=len,reverse=True)
    return res