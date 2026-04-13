#对输出序列优化

# utilities.py 正确版本（统一4空格缩进）
def colored(seq):
    # 定义碱基对应的终端颜色编码（ANSI转义码）
    bcolors = {
        'A': '\033[92m',  # 绿色
        'C': '\033[94m',  # 蓝色
        'G': '\033[93m',  # 黄色
        'T': '\033[91m',  # 红色
        'U': '\033[91m',  # 红色（RNA的U和DNA的T同色）
        'reset': '\033[0;0m'  # 重置颜色
    }
    tmpStr = ""
    for nuc in seq:
        if nuc in bcolors:
            tmpStr += bcolors[nuc] + nuc  # 给合法碱基加对应颜色
        else:
            tmpStr += bcolors['reset'] + nuc  # 非法字符重置颜色
    return tmpStr + bcolors['reset']  # 最后重置颜色，避免影响后续输出