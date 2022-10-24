def Tarantula(k :"list"):
    try:
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return (ef/(ef+nf))/(ef/(ef+nf)+ep/(ep+np))
    except:
        return -1

def sbi(k :"list"):
    try:
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return 1 - ep/(ep+ef)
    except:
        return 0
    
def Ochiai(k :"list"):
    try:
        import math
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return ef/math.sqrt((ef+ep)*(ef+nf))
    except:
        return -1
    
def Jaccard(k :"list"):
    try:
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return ef/(ef+np+nf)
    except:
        return -1
    
def ochiai2(k :"list"):
    try:
        import math
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return ef*np/math.sqrt((ef+ep)*(nf+np)*(ef+np)*(nf+ep))
    except:
        return -1

def kulc(k :"list"):
    try:
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return ef/(nf+ep)
    except:
        return -1
    
def Op2(k :"list"):
    try:
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return ef-ep/(1+ep+np)
    except:
        return -1

def Dstar2(k :"list"):
    try:
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return ef*ef/(ep+nf)
    except:
        return 0

def Dstar(k :"list"):
    try:
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return ef*ef/(ep+nf)
    except:
        return -1

def GP13(k :"list"):
    try:
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return ef*(1+1/(2*ep+ef))
    except:
        return -1

def Naish1(k :"list"):
    try:
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]

        if ef < ef + nf:
            return -1
        else:
            return ep + np - ep
    except:
        return -1

def Naish2(k :"list"):
    try:
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return ef-ep/(ep+np+1)
    except:
        return -1

def Barinel(k :"list"):
    try:
        n = k[0] + k[1] + k[2] + k[3]
        ef = k[0]
        ep = k[1]
        nf = k[2]
        np = k[3]
        return ef/(ef + ep)
    except:
        return -1


if __name__ == '__main__':
    t = [1,0,1,3]
    print(Ochiai(t))
    print(GP13(t))