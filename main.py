#DNA工具测试
from DNAToolkit import *
import random

#输入一个随机序列测试
rndDNAStr=''.join([random.choice(Nucleotides)
                   for nuc in range(50)])

print(validateSeq(rndDNAStr))
print(countFrequency(rndDNAStr))