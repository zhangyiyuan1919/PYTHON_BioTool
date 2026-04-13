#DNA工具测试
from DNAToolkit import *
from utilities import colored
import random

#输入一个随机序列测试
rndDNAStr=''.join([random.choice(Nucleotides)
                   for nuc in range(50)])

DNAStr= validateSeq(rndDNAStr)

print(f'\nSequence:{colored(DNAStr)}\n')
print(f'[1]+Sequence Length:{len(DNAStr)}\n')
print(colored(f'[2]+Nucleotide Frequency:{countFrequency(DNAStr)}\n'))
print(f'[3]+DNA/RNA Transcription:{colored(transcription(DNAStr))}\n')
print(f'[4]+Reverse Complememt:{colored(reverse_complement(DNAStr))}\n')

print(f"[5]+DNA String + Reverse Complement:\n5'{colored(DNAStr)}3'")
print(f"  {''.join(['|' for c in range(len(DNAStr))])}")
print(f"3'{colored(reverse_complement(DNAStr))}5'\n")