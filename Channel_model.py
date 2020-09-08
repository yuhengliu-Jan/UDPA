
import math
import random
import inspect
import Judgelink as JL
import numpy as np
from sympy import *


# 2019.12.3

def PL_free(Pt, fc,dist,Gt,Gr):
    # fc: 载波频率[Hz]
    # dist: 无人机之间的距离[m]
    # Gt: 发射机天线增益
    # Gr: 接收机天线增益 :
    # PL     : 路径损耗[dB]
    # nargin = inspect.getargspec(PL_free)
    fc = fc*1000000
    dist =dist*1000
    if (Gr == 0):
        nargin = 3
    else:
        nargin = 4
    lamda = 3e8/fc
    tmp = lamda/(4*math.pi*dist)
    if (nargin > 2):
        tmp = tmp * math.sqrt(Gt)
    elif (nargin > 3):
        tmp = tmp * math.sqrt(Gr)
    PL = -20 * math.log10(tmp)
    #Pt = 10*log

    Pr = Pt*pow(tmp,2)
    print('PR',Pr)
    print(10*math.log10(Pt/Pr))
    return PL
'''a1 = PL_free(1,2,0.5,1,1)
print('aaaaa',a1)'''

def PL_AG(type,fc,dn,Pt,los,nlos):
    # 基于视距非视距的路径消耗来计算无人机到地面
    # fc: 载波频率[Hz]
    # dn: 无人机与地面节点之间的距离[m]
    # type: #视距/非视距
    # los = 1
    # nlos = 2
    if (type == 0):
        PL =  20*math.log10(dn) +  20*math.log10(fc) + los
    else:
        PL = 20*math.log10(dn) +  20*math.log10(fc) + nlos
    Pr = Pt - PL

    return Pr,PL

# 基于仰角
def PL_A2G(type, theta, dist, fc):
    # fc:GHZ
    # dist : km
    if (fc == 2):
        if (theta<10):
            if (type == 0):
                PL = 20*math.log10(dist) + 20*math.log10(fc) + 92.4 + np.random.normal(0,4)
            else:
                PL = 20 * math.log10(dist) + 20 * math.log10(fc) + 92.4 + np.random.normal(
                    ((2.55 + theta) / (0.0594 + 0.0406 * theta)),
                    (-12.96 + theta) / (-1.076 + 0.078 * theta)) + np.random.normal(0, 10)
        else:
            if (type == 0):
                PL = 20 * math.log10(dist) + 20 * math.log10(fc) + 92.4 + np.random.normal(0, 4)
            else:
                PL = 20 * math.log10(dist) + 20 * math.log10(fc) + 92.4 + np.random.normal(
                    (-94.2 + theta) / (-3.44 + 0.0318 * theta), (-89.55 + theta) / (-8.87 + 0.0927 * theta)) + np.random.normal(0,10)

    elif (fc == 3.5):
        if (theta<10):
            if (type == 0):
                PL = 20*math.log10(dist) + 20*math.log10(fc) + 92.4 + np.random.normal(0,4)
            else:
                PL = 20 * math.log10(dist) + 20 * math.log10(fc) + 92.4 + np.random.normal(((2.7 +theta)/(0.059+0.0376*theta)),(-12.24+theta)/(-1.006+0.0788*theta)) + np.random.normal(0,10)
        else:
            if (type == 0):
                PL = 20 * math.log10(dist) + 20 * math.log10(fc) + 92.4 + np.random.normal(0, 4)
            else:
                PL = 20 * math.log10(dist) + 20 * math.log10(fc) + 92.4 + np.random.normal(
                    ((-92.9 + theta) / (-3.14 + 0.0302 * theta)), (-89.06 + theta) / (-8.63 + 0.0921 * theta)) + np.random.normal(0,10)

    elif (fc == 5.5):
        if (theta<10):
            if (type == 0):
                PL = 20*math.log10(dist) + 20*math.log10(fc) + 92.4 + np.random.normal(0,4)
            else:
                PL = 20 * math.log10(dist) + 20 * math.log10(fc) + 92.4 + np.random.normal(((2.636+theta)/(0.0554+0.0352*theta)),(-12.4+theta)/(-0.998+0.0769*theta)) + np.random.normal(0,10)
        else:
            if (type == 0):
                PL = 20 * math.log10(dist) + 20 * math.log10(fc) + 92.4 + np.random.normal(0, 4)
            else:
                PL = 20 * math.log10(dist) + 20 * math.log10(fc) + 92.4 + np.random.normal(
                    ((-92.8 + theta) / (-2.955 + 0.0285 * theta)), (-89.54 + theta) / (-8.474 + 0.09 * theta)) + np.random.normal(0,10)

    return PL





def PL_exAG(Pt, dist , Au , type):
    # 基于路径衰减指数是根据PPP和仰角的无人机到地面
    # Au = 2

    aaf = 20
    if (type == 0) :
        Pr = Pt * dist^(Au)
    else:
        Pr = aaf * Pt * dist^(Au)

    return Pr


def PL_G2A(Pt, fc, d, d0, n, sigma):
    # d : 基站和移动台之间的距离[m]
    # d0 : 参考距离[m] 100
    # n : 路径损耗指数    # sigma : 方差[dB]
    lamda = 3e8 / fc
    PL = -20 * math.log10(lamda/(4*math.pi*d0)) + 10 * n *math.log10(d/d0)
    if (sigma == 0):
        Pr = Pt - PL
    else:
        PL = PL + sigma * random.uniform(0,1)
        Pr = math.log10(Pt) - math.log10(pow(10,PL))

    return Pr

def Pr_G2G(Pd, dist, type, g0):
    #type 环境
    if (type == 1):
        ad = -3
    elif(type == 2):
        ad = -4
    elif (type == 3):
        ad = -6
    Prd = Pd * pow(dist,ad) * g0
    return  Prd

# 2019.12.17
def cover_G2G(pt, fc,type,g0,minpr, eps):
    print('地面节点的通信范围')
    if (type == 1):
        ad = -3
    elif(type == 2):
        ad = -4
    elif (type == 3):
        ad = -6

    dis = pow(pt/(minpr*g0),ad)
    print(pt * pow(20,ad) * g0)
    return dis

def log_G2G(fc,d,d0,n,Xa):
    # xa = 3
    lamda = 3e8 / fc
    pl = -20*math.log10(lamda/(4*math.pi*d0))+10*n*math.log10(d/d0) + Xa*np.random.random()*d
    print('pl',pl)
    return pl

'''ss = log_G2G(2000000,50,1000,3,3)

xx = cover_G2G(1,2,1,1,1e-7,-7)
print('cover_g2g',xx)'''




# 重新写地面节点之间的接收信号强度


def G2G_SINR(Pd):
    return 0


def panduan(pr,eps):
    temp = []
    for i in pr:
        if (i>=eps):
            temp.append(0)
        else:
            temp.append(1)
    return temp

#2019.12.19 有点问题
def tan_n0(tran_node, Pt, ad, gi):
    Id = 0
    temp = []
    for i in range(tran_node.shape[0]):
        for j in range(tran_node.shape[0]):
            if (i == j):
                print('是自己人')
            else:
                n1 = JL.point(tran_node[i][0],tran_node[i][1],tran_node[i][2])
                # print(n1)
                n2 = JL.point(tran_node[j][0],tran_node[j][1],tran_node[j][2])
                di = JL.two_node_disl(n1,n2)
                Id = Pt*pow(di,ad)*gi[j]+Id
        temp.append(Id)
        Id = 0
    return  temp

# 别的节点干扰
def oth_n(tran_node, sum_list):
    print('我还没写')


# 还需要补全改进
def cover_dis(eps, Pt,  los, nlos,type, theta, dist, fc):
    print('计算距离')
    dist = 0.2
    Ngw = -174
    N0 = 10
    # los = 2
    # nlos = 36.54
    if (type == 0):
        PL = 20 * math.log10(dist) + 20 * math.log10(fc) + 92.4 + los
    else:
        PL = 20 * math.log10(dist) + 20 * math.log10(fc) + 92.4 + nlos
    print(PL)
    p = 95
    sn = pow(10,(p - 20*math.log10(fc)-92.4+1)/20)
    print('sn',sn)
    print(math.cos((42/180)*math.pi))
    dn = math.sqrt(pow(sn,2) * pow(math.cos(((42/180)*math.pi)),2))
    print('dn',dn)
    hn =  math.sqrt(pow(sn,2) * pow(math.sin(((42/180)*math.pi)),2))
    print('hn',hn)
    # Pr = Pt/pow(10,PL/10)

    x = Symbol('x')
    # cdd = solve(((10*log(Pt*1000))/pow(10,(20 * log(x) + 20 * log(fc) + 92.4 + los)/10)) +70 ,x)
    max = 10*math.log10((Pt*1000)/1e-7)
    print(solve((20 * log(x) + 20 * log(fc) + 92.4 + los)- max,x))
    print(solve((20 * log(x) + 20 * log(fc) + 92.4 + nlos) - max, x))
    return max

    #print('cdd',cdd*1000)

def get_build(mobile_file_path):
    build_list = []
    with open(mobile_file_path, 'r') as f:
        for line in f :
            line_list = re.split('[\s]', line)
            # build_list.append((int(float(line_list[0])), int(float(line_list[1]))))
            build_list.append((float(line_list[0]), float(line_list[1]), float(line_list[0]), float(line_list[4]), float(line_list[4]), float(line_list[6]), float(line_list[0]), float(line_list[6]), float(line_list[7])))
            # print(line_list[7])
        # print(build_list)
    return build_list


