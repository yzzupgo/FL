# totalfailed indicates the number of failed test cases for the original program P
# type = int,tf 在原程序上faile的数量
# type = list,f p 在变异体mij上f2p 和p2f的数量
# type = int ,f2p p2f 所有变异体上f2p 和p2f的数量
import math
type_mbfl = 'max'
alpha = 0.5

def select_type(sus_list):
    sus_list_format = list(map(lambda x: float('%.5f' % x),sus_list))
    if type_mbfl == 'max':
        return max(sus_list_format)
    elif type_mbfl == 'ave':
        return sum(sus_list_format)/len(sus_list_format)
    elif type_mbfl == 'none':
        frequency_dict = dict()
        for sus in sus_list_format:
            if sus not in frequency_dict:
                frequency_dict[sus] = 0
            frequency_dict[sus] += 1
        return max(sorted(list(map(lambda x: [frequency_dict[x], x],
                                   list(frequency_dict.keys()))),
                          key=lambda x: x[1], reverse=True))[1]
    elif type_mbfl == 'frequency':
        return sorted(sus_list_format, reverse=True)


def nsus(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    akf_list, akp_list, anf_list, anp_list = touple[6], touple[7], touple[8], touple[9]
    sus = []
    for i in range(len(f)):
        try:
            sus.append(f[i]/tf-alpha*p[i]/tp)
        except:
            sus.append(-1)
    return select_type(sus)

def muse(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    akf_list, akp_list, anf_list, anp_list = touple[6], touple[7], touple[8], touple[9]
    sus = []
    for i in range(len(f)):
        try:
            sus.append((akf_list[i]+anf_list[i]) - f2p/p2f*(akp_list[i]+anp_list[i]))
        except:
            sus.append(-1)
    return select_type(sus)


def metallaxis(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    sus = []
    for i in range(len(f)):
        try:
            sus.append(f[i]/math.sqrt(tf*(f[i]+p[i])))
        except:
            sus.append(-1)
    return select_type(sus)


def Tarantula(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    akf_list, akp_list, anf_list, anp_list = touple[6], touple[7], touple[8], touple[9]

    sus = []
    for i in range(len(akf_list)):
        akf = akf_list[i]
        akp = akp_list[i]
        anf = anf_list[i]
        anp = anp_list[i]
        try:
            sus.append((akf/(akf+anf))/((akf/(akf+anf))+(akp/(akp+anp))))
        except:
            sus.append(-1)

    return select_type(sus)

def Op2(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    akf_list, akp_list, anf_list, anp_list = touple[6], touple[7], touple[8], touple[9]

    sus = []
    for i in range(len(akf_list)):
        akf = akf_list[i]
        akp = akp_list[i]
        anf = anf_list[i]
        anp = anp_list[i]
        try:
            sus.append(akf-akp/(akp+anp+1))
        except:
            sus.append(-1)

    return select_type(sus)

def Ochiai(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    akf_list, akp_list, anf_list, anp_list = touple[6], touple[7], touple[8], touple[9]

    sus = []
    for i in range(len(akf_list)):
        akf = akf_list[i]
        akp = akp_list[i]
        anf = anf_list[i]
        anp = anp_list[i]
        try:
            sus.append(akf/math.sqrt((akf+anf)*(akf+akp)))
        except:
            sus.append(-1)

    return select_type(sus)

def Dstar(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    akf_list, akp_list, anf_list, anp_list = touple[6], touple[7], touple[8], touple[9]

    sus = []
    for i in range(len(akf_list)):
        akf = akf_list[i]
        akp = akp_list[i]
        anf = anf_list[i]
        anp = anp_list[i]
        try:
            sus.append(akf*akf/(akp+anf))
        except:
            sus.append(-1)

    return select_type(sus)

def GP13(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    akf_list, akp_list, anf_list, anp_list = touple[6], touple[7], touple[8], touple[9]

    sus = []
    for i in range(len(akf_list)):
        akf = akf_list[i]
        akp = akp_list[i]
        anf = anf_list[i]
        anp = anp_list[i]
        try:
            sus.append(akf*(1+1/(2*akp+akf)))
        except:
            sus.append(-1)

    return select_type(sus)

def Naish1(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    akf_list, akp_list, anf_list, anp_list = touple[6], touple[7], touple[8], touple[9]

    sus = []
    for i in range(len(akf_list)):
        akf = akf_list[i]
        akp = akp_list[i]
        anf = anf_list[i]
        anp = anp_list[i]
        if akf < tf:
            sus.append(-1)
        else:
            sus.append(tp-akp)

    return select_type(sus)

def Naish2(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    akf_list, akp_list, anf_list, anp_list = touple[6], touple[7], touple[8], touple[9]

    sus = []
    for i in range(len(akf_list)):
        akf = akf_list[i]
        akp = akp_list[i]
        anf = anf_list[i]
        anp = anp_list[i]
        sus.append(akf-akp/(akp+anp+1))

    return select_type(sus)

def Jaccard(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    akf_list, akp_list, anf_list, anp_list = touple[6], touple[7], touple[8], touple[9]

    sus = []
    for i in range(len(akf_list)):
        akf = akf_list[i]
        akp = akp_list[i]
        anf = anf_list[i]
        anp = anp_list[i]
        try:
            sus.append(akf/(akf+anf+anp))
        except:
            sus.append(-1)

    return select_type(sus)

def Barinel(touple):
    tf, tp, f, p, f2p, p2f = touple[0], touple[1], touple[2], touple[3], touple[4], touple[5]
    akf_list, akp_list, anf_list, anp_list = touple[6], touple[7], touple[8], touple[9]

    sus = []
    for i in range(len(akf_list)):
        akf = akf_list[i]
        akp = akp_list[i]
        anf = anf_list[i]
        anp = anp_list[i]
        try:
            sus.append(akf/(akf+akp))
        except:
            sus.append(-1)

    return select_type(sus)


# anf 存活 执行失败
# anp 存活 执行通过
# akf 杀死 失败 f2p
# akp 杀死 通过 p2f
#
# 计算单个变异体怀疑度。
# 通过单个变异体怀疑度取平均获取语句怀疑度。