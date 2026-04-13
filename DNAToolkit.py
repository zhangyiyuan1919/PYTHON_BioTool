#DNA工具包

import collections
from structures import *

# 工具1：定义合法DNA碱基（支持IUPAC简并碱基可自行扩展）
def validateSeq(dna_seq):
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

#工具2：统计碱基数量
def countFrequency(seq):
    tmpFreqDict={"A":0,"T":0,"C":0,"G":0} #计数字典
    for nuc in seq:
        tmpFreqDict[nuc]+=1
    return tmpFreqDict
#另一种解法
#def countNucFrequency(seq):
#    return dict(collections.Counter(seq))

#工具3：序列转录 DNA→RNA
def transcription(seq):
    """DNA→RNA转录，用U代替T"""  #在自定工具包中相当于加入注释
    return seq.replace("T","U")

#工具4：互补序列
def reverse_complement(seq):
    """DNA的互补链，顺序是5‘→3’"""
    return ''.join([DNA_ReverseComplement[Nuc] for Nuc in seq])[::-1]#[::-1]首位调转互补后的序列