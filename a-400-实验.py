import time
import Get_Move as Gm
import numpy as np
import Channel_model as CM
import math
import Judgelink as JL
import networkx as nx
import copy
import hull
import alfour
from scipy import stats
from pyscipopt import Model, quicksum, multidict, SCIP_PARAMSETTING, exp, log, sqrt
import matplotlib.pyplot as plt
from sympy import *


def broadcast1(uav, node_list, type, eps, fc, Pt):
    node_num = []
    node_inti = []
    node_inf = []
    mmm = 0
    mum = 0
    qoss = 0
    shiji = []
    for node in node_list:
        temp = 0
        k = 1
        s = 0
        aim = JL.point(uav[0], uav[1], uav[2])
        BB = JL.point(node[0, 1], node[0, 2], node[0, 3])
        dna = pow((aim.x - BB.x), 2) + pow((aim.y - BB.y), 2) + pow((aim.z - BB.z), 2)
        dn = math.sqrt(dna)
        temp_jiao = uav[2] / dn
        if (temp_jiao > 1):
            temp_jiao = 0.99
        theta = math.degrees(math.asin(temp_jiao))
        A = JL.J_broadcast(aim, BB)

        if (A == 1):
            shiji.append((node[0, 0], node[0, 1], node[0, 2], dn))
            if (dn < 1068):
                qoss = qoss + 1
                node_inf.append((node[0], theta, dn, 1, uav))
                node_num.append(node[0, 0])
                node_inti.append((node[0, 1], node[0, 2]))



    return node_inti, node_inf, node_num, shiji

def bestpp(uav, node_list):
    node_num = []
    tnn = []
    node_inti = []
    node_inf = []
    qoss = 0
    shiji = []
    alloss = []
    Pt = 1000
    eps = -70
    fc = 2
    for node in node_list:
        aim = JL.point(uav[0], uav[1], uav[2])
        BB = JL.point(node[0, 1], node[0, 2], node[0, 3])
        dna = pow((aim.x - BB.x), 2) + pow((aim.y - BB.y), 2) + pow((aim.z - BB.z), 2)
        dn = math.sqrt(dna)
        temp_jiao = uav[2] / dn
        if (temp_jiao > 1):
            temp_jiao = 0.99
        theta = math.degrees(math.asin(temp_jiao))
        A = JL.J_broadcast(aim, BB)

        if (A == 1):
            if (dn < 1068):
                qoss = qoss + 1
                shiji.append(node[0, 0])
                node_inf.append((node[0], theta, dn, 1, uav))
                node_inti.append((node[0, 1], node[0, 2]))
                tnn.append(node[0, 0])

            alloss.append(node[0, 0])

            PLc = CM.PL_A2G(0, theta, dn / 1000, fc)
            PL = 10 * math.log10(Pt) - PLc
            if (PL > eps):
                tnn.append(node[0, 0])

    node_num = list(set(tnn))

    return node_inti, node_inf, node_num, shiji, alloss


def pre_grid(D):
    bian = 20
    tnum = int((D * 1000) / bian)
    tt = tnum * tnum
    s = 1
    ditu = []
    for kk in range(tt):
        ditu.append(int(s))
        s = s + 1
    gird_list = []
    org_gird = np.array([bian / 2, bian / 2])
    add_r = np.array([0, bian])
    add_l = np.array([bian, 0])
    num = tnum - 1
    for i in range(tnum):
        if (i == 0):
            temp = org_gird
            temp1 = temp
            gird_list.append(temp)
        else:
            temp = temp + add_l
            temp1 = temp
            gird_list.append(temp)
            # print(temp)
        for j in range(num):
            temp1 = temp1 + add_r
            gird_list.append(temp1)
    fina = []
    for xxe in gird_list:
        fina.append((xxe[0], xxe[1], 1))

    return bian, tnum, ditu, fina


def is_in_line(point, o, d):
    # 先判断该点是否在线段范围内,如果不在,
    # 则就算在该方程的直线上但也不在该线段上
    if o[1] > point[1] and d[1] > point[1]:
        return False
    if o[1] < point[1] and d[1] < point[1]:
        return False
    if o[0] > point[0] and d[0] > point[0]:
        return False
    if o[0] < point[0] and d[0] < point[0]:
        return False

    if o[0] == d[0]:
        return True if point[0] == o[0] else False

    # 通过输入的两个点计算一元一次方程,通过输入x,计算y
    a = (d[1] - o[1]) / (d[0] - o[0])
    b = (d[0] * o[1] - o[0] * d[1]) / (d[0] - o[0])
    y = a * point[0] + b
    return True if y == point[1] else False


def is_ray_intersects_segment(point, o, d):
    if o[1] == d[1]:  # 如果线段是水平直线,则直接无交点
        return False
    # 先判断该点是否在线段范围内,如果不在,
    # 则就算在该方程的直线上但也不在该线段上
    if o[1] > point[1] and d[1] > point[1]:  # 该点纵坐标小于等于线段最小纵坐标
        return False
    if o[1] < point[1] and d[1] < point[1]:  # 该点纵坐标大于等于线段最大纵坐标
        return False

    # 利用三角形比例关系求交点
    x = d[0] - (d[0] - o[0]) * (d[1] - point[1]) / (d[1] - o[1])

    # 如果该点的横坐标小于相同水平线上交点的横坐标,则有交点
    return True if point[0] < x else False


def is_point_in_polygon(point_coord, polygon_coords, is_contains_edge):
    intersect_count = 0  # 交点个数
    for polygon in polygon_coords:
        # 循环每条边
        for i in range(len(polygon) - 1):
            origin_point = polygon[i]
            destination_point = polygon[i + 1]

            # 是否包含存在直线上的点
            if is_in_line(point_coord, origin_point, destination_point):
                return True if is_contains_edge else False

            if is_ray_intersects_segment(point_coord, origin_point, destination_point):
                # 有交点就加1
                intersect_count += 1

    # 如果恰好与端点相交,则查看相应的端点,并减去与该点相同纵坐标且在该点右侧的点的个数
    endpoint_intersects_count = 0
    for polygon in polygon_coords:
        # 遍历每一个点,但是因为最后一个点和第一个点是重复的,所以最后一个不遍历
        for i in polygon[:-1]:
            if i[1] == point_coord[1] and i[0] > point_coord[0]:
                endpoint_intersects_count += 1
    intersect_count -= endpoint_intersects_count
    return True if intersect_count % 2 == 1 else False


def uav_direction_ang(uav, inf):
    uav_o = np.array([uav[0], uav[1]])
    origin = np.array([uav[0] + 1, uav[1]])
    node_pos = []
    si = 0
    co = 0
    for c_node in inf:
        # node_pos = np.array([node_pos[1],node_pos[2]])
        temp = np.array([c_node[1], c_node[2]])
        x = origin - uav_o
        y = temp - uav_o
        Lx = np.sqrt(x.dot(x))
        Ly = np.sqrt(y.dot(y))
        cos_angle = x.dot(y) / (Lx * Ly)
        angle = np.arccos(cos_angle)
        # print(angle)
        if (y[0] >= 0 and y[1] >= 0):
            angle1 = angle
        elif (y[0] < 0 and y[1] >= 0):
            angle1 = angle
        else:
            angle1 = (2 * math.pi) - angle
        node_pos.append(angle1)
        si = si + math.sin(angle1)
        co = co + math.cos(angle1)
    # print('xxxx',co/si)
    return node_pos


def girl_ang(uav, inf):
    uav_o = np.array([uav[0], uav[1]])
    origin = np.array([uav[0] + 1, uav[1]])
    node_pos = []
    for c_node in inf:
        # node_pos = np.array([node_pos[1],node_pos[2]])
        temp = np.array([c_node[0], c_node[1]])
        x = origin - uav_o
        y = temp - uav_o
        Lx = np.sqrt(x.dot(x))
        Ly = np.sqrt(y.dot(y))
        cos_angle = x.dot(y) / (Lx * Ly)
        angle = np.arccos(cos_angle)
        # print(angle)
        if (y[0] >= 0 and y[1] >= 0):
            angle1 = angle
        elif (y[0] < 0 and y[1] >= 0):
            angle1 = angle
        else:
            angle1 = (2 * math.pi) - angle
        node_pos.append((c_node[0], c_node[1], angle1))
    return node_pos


def compute_Attract11(X, Xsum, k, neib):
    # 输入参数为当前坐标，目标坐标，增益常数,分量和力的角度
    # 把路径上的临时点作为每个时刻的Xgoal
    uav_o = np.array([X[0], X[1]])
    Yatx = 0
    Yaty = 0
    for node in Xsum:
        R = pow((X[0] - node[0]), 2) + pow((X[1] - node[1]), 2)  # + pow(X[3] - Xsum[0][3], 2) #%路径点和目标的距离平方
        r = math.sqrt(R)  # 路径点和目标的距离
        Yatx = k * r * math.cos(node[2]) + Yatx  # angle=Y(1)
        Yaty = k * r * math.sin(node[2]) + Yaty
    uuu = math.sqrt(Yaty * Yaty)
    dire = math.asin(uuu / math.sqrt(pow(Yatx, 2) + pow(Yaty, 2)))
    if (Yatx > 0 and Yaty > 0):
        dire = dire
    elif (Yatx < 0 and Yaty > 0):
        dire = math.pi - dire
    elif (Yatx < 0 and Yaty < 0):
        dire = math.pi + dire
    else:
        dire = math.pi + math.pi - dire
    tx1 = math.cos(dire)
    ty1 = math.sin(dire)
    XX = (tx1 * 0.7)
    YY = (ty1 * 0.7)

    neix = 0
    neiy = 0
    if (neib != []):
        for tt in neib:
            temp = np.array([tt[0], tt[1]])
            x = np.array([1, 0])
            y = temp - uav_o
            Lx = np.sqrt(x.dot(x))
            Ly = np.sqrt(y.dot(y))
            cos_angle = x.dot(y) / (Lx * Ly)
            angle = np.arccos(cos_angle)
            # print(angle)
            if (y[0] >= 0 and y[1] >= 0):
                angle1 = angle
            elif (y[0] < 0 and y[1] >= 0):
                angle1 = angle
            else:
                angle1 = (2 * math.pi) - angle
            neix = k * Ly * math.cos(angle1) + neix
            neiy = k * Ly * math.sin(angle1) + neiy
        uuu = math.sqrt(neiy * neiy)
        dire = math.asin(uuu / math.sqrt(pow(neix, 2) + pow(neiy, 2)))
        if (neix > 0 and neiy > 0):
            dire = dire
        elif (neix < 0 and neiy > 0):
            dire = math.pi - dire
        elif (neix < 0 and neiy < 0):
            dire = math.pi + dire
        else:
            dire = math.pi + math.pi - dire
        nx1 = math.cos(dire)
        ny1 = math.sin(dire)
        XX = (tx1 * 0.7) + (nx1 * 0.3)
        YY = (ty1 * 0.7) + (ny1 * 0.3)

    uuu = math.sqrt(pow(YY, 2))
    dire = math.asin(uuu / math.sqrt(pow(XX, 2) + pow(YY, 2)))
    if (XX > 0 and YY > 0):
        dire = dire
    elif (XX < 0 and YY > 0):
        dire = math.pi - dire
    elif (XX < 0 and YY < 0):
        dire = math.pi + dire
    else:
        dire = math.pi + math.pi - dire

    return dire


def uav_deal_ang(uav, inf):
    uav_o = np.array([uav[0], uav[1]])
    node_pos = []
    for c_node in inf:
        # node_pos = np.array([node_pos[1],node_pos[2]])
        temp = np.array([c_node[0, 1], c_node[0, 2]])
        origin = np.array([temp[0] + 1, temp[1]])
        x = origin - temp
        y = uav_o - temp
        Lx = np.sqrt(x.dot(x))
        Ly = np.sqrt(y.dot(y))
        cos_angle = x.dot(y) / (Lx * Ly)
        angle = np.arccos(cos_angle)
        # print(angle)
        if (y[0] >= 0 and y[1] >= 0):
            angle1 = angle
        elif (y[0] < 0 and y[1] >= 0):
            angle1 = angle
        else:
            angle1 = (2 * math.pi) - angle
        node_pos.append((c_node[0, 0], angle1, c_node[0, 1], c_node[0, 2]))

    return node_pos


def kenengxing2(dist):
    fc = 2
    theta = math.degrees(math.atan(0.4 / dist))
    PL = (20 * math.log10(math.sqrt(0.16 + pow(dist, 2)))) + (20 * math.log10(fc)) + 92.4
    ca = 100 - PL
    if (theta > 10):
        g = -94.2
        h = -3.44
        i = 0.0318
        g1 = -89.55
        h1 = -8.87
        i1 = 0.0927
    else:
        g = 2.55
        h = 0.0594
        i = 0.0406
        g1 = -12.96
        h1 = -1.076
        i1 = 0.078

    mu = (g + theta) / (h + (i * theta))
    sigma = (g1 + theta) / (h1 + (i1 * theta))
    sigma = math.sqrt(pow(sigma, 2) + 100)

    Plos = (120 - (120 / (1 + pow(theta / 24.3, 1.229)))) / 100

    glos = 1-stats.norm.cdf(-ca, 0, 4)
    nglos = 1 - glos

    final = Plos * nglos + 1 - Plos

    return final


def UAV_Information_collection(drone, node_list, girl_list, neib, enenti, ntemp, aall, ini):
    tempn1 = []
    tempna = []
    node_inf1 = []
    for tp1 in ntemp:
        tempn1.append(copy.deepcopy(tp1))
    for tp2 in aall:
        tempna.append(copy.deepcopy(tp2))
    for tp3 in ini:
        node_inf1.append(copy.deepcopy(tp3))

    attract = []
    ag = 0
    while 1:
        ag1 = []
        ag2 = []
        ag3 = []
        ag4 = []
        for ainode in girl_list:
            if (ainode[0] < 500 and ainode[1] > 500):
                ag1.append(ainode)
            if (ainode[0] < 500 and ainode[1] < 500):
                ag2.append(ainode)
            if (ainode[0] > 500 and ainode[1] < 500):
                ag3.append(ainode)
            if (ainode[0] > 500 and ainode[1] > 500):
                ag4.append(ainode)

        if (ag1 != []):
            attract = ag1
            break
        if (ag2 != [] and ag1 == []):
            attract = ag2
            break
        if (ag3 != [] and ag2 == [] and ag1 == []):
            attract = ag3
            break
        if (ag4 != [] and ag2 == [] and ag1 == [] and ag3 == []):
            attract = ag4
            break
        if (ag2 == [] and ag1 == [] and ag3 == [] and ag4 == []):
            ag = 1
            break
    if (ag == 0):
        if (attract == []):
            print('xxxxx')

        girlist = girl_ang(drone, attract)
        FFF = compute_Attract11(drone, girlist, 1, neib)
        # FFF = compute_Attract12(drone, girlist, 1)
        addx = 161 * math.cos(FFF)
        addy = 161 * math.sin(FFF)
        drone = [drone[0] + addx, drone[1] + addy, 400]
        print('位置', drone)

        node_inti, node_inf, node_num, shiji = broadcast1(drone, node_list, 0, -70, 2, 1000)
        print('...........................分割线........................')
        node_inf1.append(copy.deepcopy(node_inf))
        # print('数据',enenti,len(node_inf1),node_inf1)
        tempn1.append((drone, len(node_num), copy.deepcopy(node_num)))

        k = 0
        for dtemp in girl_list:
            if (dtemp[0] != 500 and dtemp[1] != 500):
                tedis = math.sqrt(pow(dtemp[0] - drone[0], 2) + pow(dtemp[1] - drone[1], 2)) / 1000
                if (tedis < 0.4):
                    dddd = kenengxing2(tedis)
                    temp = dtemp[2] * dddd
                    girl_list[k] = [dtemp[0], dtemp[1], temp]
                    if (girl_list[k][2] < 0.02):
                        girl_list[k] = [500, 500]
            k = k + 1

        # print('记录', enenti, FFF, drone)
        enenti = enenti + 1
        UAV_Information_collection(drone, node_list, girl_list, neib, enenti, tempn1, tempna, node_inf1)
    if (ag == 1):
        sx = node_inf1[2][0][4][0] - drone[0]
        sy = node_inf1[2][0][4][1] - drone[1]
        stp = math.sqrt(pow(sx, 2) + pow(sy, 2)) / 161
        sx = sx / stp
        sy = sy / stp
        for tea in range(int(stp)):
            drone = [drone[0] + sx, drone[1] + sy, 300]
            print('位置', drone)
            node_inti, node_inf, node_num, shiji = broadcast1(drone, node_list, 0, -70, 2, 1000)
            print('...........................分割线........................')
            node_inf1.append(copy.deepcopy(node_inf))
            # print('数据',enenti,len(node_inf1),node_inf1)
            tempn1.append((drone,len(node_num), copy.deepcopy(node_num)))
            enenti = enenti + 1
        jjlu = []
        ccnode = []
        for ii in tempn1:
            for jj in ii[2]:
                if (jj not in ccnode):
                    ccnode.append(jj)

        ccnode.sort()
        print('ccccc', len(ccnode), ccnode)
        for shengyunn in range(len(node_list)):
            if (shengyunn not in ccnode):
                print('没扫描到：', shengyunn)

        deal_all_inf(node_inf1, node_list, ccnode,tempn1)


def deal_all_inf(node_inf, init_position_matrix, alnodes,hisl):
    # 将邻居，同一方向的夹角考虑进去，还有无人机的步长。2020.3.2

    time = len(node_inf)
    print('time', time)
    print(node_inf[0])
    uav_pos = []
    node_inf1 = []
    for ii in range(time):
        if(node_inf[ii] == []):
            uav_pos.append(hisl[ii][0])
        else:
            uav_pos.append(node_inf[ii][0][4])
    # print(uav_pos)

    for evet in range(time):
        node_time = []
        node_idd = []
        node_los = []
        for ntt in node_inf[evet]:
            if (ntt[0][0, 0] not in node_idd):
                node_idd.append(ntt[0][0, 0])
                node_time.append(ntt[0][0])
                node_los.append(ntt[3])
        node_ang = uav_deal_ang(uav_pos[evet], node_time)
        ttemp = []
        for ttt in range(len(node_los)):
            ttemp.append((node_ang[ttt][0], node_ang[ttt][1], node_los[ttt], node_ang[ttt][2], node_ang[ttt][3]))
        node_inf1.append(ttemp)
    print(node_inf1[0])

    # 角度范围
    angcc = {i: [] for i in range(len(init_position_matrix))}
    allang = {j: [] for j in range(len(init_position_matrix))}
    alang = {k: [] for k in range(len(init_position_matrix))}
    for en in range(len(init_position_matrix)):
        if (en not in alnodes):
            alang[en] = []
        else:
            temp = []
            temp1 = []
            ss = 0
            for iii in node_inf1:
                m = 0
                for jjj in iii:
                    if (jjj[0] == en):
                        temp.append((jjj[1], jjj[2], jjj[3], jjj[4]))
                        temp1.append(jjj[1])
                        m = 1
                if (m == 0 and en in alnodes):
                    pdn = pow(uav_pos[ss][0] - init_position_matrix[en][0, 1], 2) + pow(
                        uav_pos[ss][1] - init_position_matrix[en][0, 2], 2)
                    if (pdn >= 980624):
                        xxx = uav_deal_ang(uav_pos[ss], [init_position_matrix[en]])
                        temp.append((xxx[0][1], -1, xxx[0][2], xxx[0][3]))
                        temp1.append(xxx[0][1])
                    else:
                        xxx = uav_deal_ang(uav_pos[ss], [init_position_matrix[en]])
                        temp.append((xxx[0][1], 0, xxx[0][2], xxx[0][3]))
                        temp1.append(xxx[0][1])
                ss = ss + 1
            angcc[en] = temp
            allang[en] = copy.deepcopy(temp1)
            temp1.sort()
            alang[en] = temp1

    ob_map = {m: [] for m in range(len(init_position_matrix))}
    for init in range(len(init_position_matrix)):

        temppx = []
        for inix in range(len(alang[init])):
            if (inix == len(alang[init]) - 1):
                tempp = (alang[init][inix] - alang[init][inix - 1]) / 2
                xx1 = alang[init][inix] - tempp
                xx3 = alang[init][inix] + tempp
                if (xx3 > (math.pi * 2)):
                    xx3 = math.pi * 2
            elif (inix == 0):
                tempp = (alang[init][1] - alang[init][0]) / 2
                xx1 = alang[init][0] - tempp
                if (alang[init][0] == 0):
                    xx1 = 0
                elif (xx1 < 0):
                    xx1 = 0
                xx3 = alang[init][0] + tempp
            else:
                tempp = (alang[init][inix] - alang[init][inix - 1]) / 2
                tempp1 = (alang[init][inix + 1] - alang[init][inix]) / 2
                xx1 = alang[init][inix] - tempp
                xx3 = alang[init][inix] + tempp1
            temppx.append((xx1, alang[init][inix], xx3))
        ob_map[init] = temppx
    # print('1212313',len(ob_map[0]),ob_map[0])

    zob_map = {z: [] for z in range(len(init_position_matrix))}
    for fin in range(len(angcc)):
        tttemp = []
        for xxf in range(len(angcc[fin])):
            jisu = 0
            for mm in angcc[fin]:
                if (ob_map[fin][xxf][1] == mm[0]):
                    diss = math.sqrt(pow(mm[2] - uav_pos[jisu][0], 2) + pow(mm[3] - uav_pos[jisu][1], 2))
                    zhi = math.atan(398.5 / diss)
                    tttemp.append((jisu, zhi, mm[1], ob_map[fin][xxf][0], ob_map[fin][xxf][1], ob_map[fin][xxf][2]))
                    break
                jisu = jisu + 1
        zob_map[fin] = tttemp
    '''
    cn1 = 0
    for ixi in node_num:
        ix = int(ixi)
        if (pow(init_position_matrix[ix][0, 1] - uav_pos[time - 1][0], 2) + pow(
                init_position_matrix[ix][0, 2] - uav_pos[time - 1][1],
                2) < 640000):
            cn1 = cn1 + 1
    print('少少少：', cn1)
    '''
    print("------------------------split line------------------------")

    demo = alfour.PSO(alnodes, init_position_matrix, zob_map, 10)
    demo.fanwei()
    demo.initial()
    demo.solving(15)
    result = demo.returnbest()
    result.sort(key=takeFirst, reverse=True)
    locm = [result[0][1], result[0][2], 400]
    print('多多少少：', len(node_num), node_num)
    print(result[0])

    adaptive_dep(uav_pos, zob_map, init_position_matrix, locm, alnodes, 10, hisl)



def adaptive_dep(uav_pos, zob_map, init_position_matrix, locm, alnodess, ctime, hisl):
    suoyou = []
    for tp2 in hisl:
        suoyou.append(copy.deepcopy(tp2))
    print('ctime', ctime)
    ctime = ctime - 1
    node_inti, node_inf2, node_num, shiji, alloss = bestpp(locm, init_position_matrix)
    suoyou.append((locm,len(node_num),node_num))
    print('范围内',len(node_num))
    print('实际符合：', len(shiji))
    # print(node_num)
    '''
    losenode = []
    for llos in shiji:
        if (llos[0] not in node_num):
            losenode.append(llos)
    print('缺少的点：', losenode)
    '''
    newnodes = []
    for kk in node_num:
        if (kk not in alnodess):
            alnodess.append(kk)
            newnodes.append(kk)
    alnodess.sort()
    jj = len(uav_pos) - 1
    # print('目前搜索到的节点数', len(alnodess))

    newml = {k: [] for k in range(len(init_position_matrix))}

    for q1 in newnodes:
        q1 = int(q1)
        tempp1 = []
        tempp2 = []
        for qqq1 in range(len(uav_pos)):
            qq1 = uav_deal_ang(uav_pos[qqq1], [init_position_matrix[q1]])
            dab = math.sqrt(pow(init_position_matrix[q1][0, 1] - uav_pos[qqq1][0], 2) + pow(
                init_position_matrix[q1][0, 2] - uav_pos[qqq1][1], 2))
            zab = math.atan(398.5 / dab)
            if (dab < 990.2):
                tempp1.append((qqq1, qq1[0][1], zab, 0))
            else:
                tempp1.append((qqq1, qq1[0][1], zab, -1))
        uav_pos.append(locm)
        tempp1.sort(key=takeSecond)
        print('ttttt', len(tempp1), tempp1)
        for inix in range(len(tempp1)):

            if (inix == len(tempp1) - 1):
                tempp = (np.array(tempp1[inix][1]) - np.array(tempp1[inix - 1][1])) / 2
                xx1 = tempp1[inix][1] - tempp
                xx3 = tempp1[inix][1] + tempp
                if (xx3 > (math.pi * 2)):
                    xx3 = math.pi * 2
            elif (inix == 0):
                tempp = (np.array(tempp1[1][1]) - np.array(tempp1[0][1])) / 2
                tt = copy.deepcopy(tempp1[0][1])
                xx1 = tempp1[0][1] - tempp
                if (tempp1[0] == 0):
                    xx1 = 0
                elif (xx1 < 0):
                    xx1 = 0
                xx3 = tempp1[0][1] + tempp
            else:
                tempp = (np.array(tempp1[inix][1]) - np.array(tempp1[inix - 1][1])) / 2
                tempp12 = (np.array(tempp1[inix + 1][1]) - np.array(tempp1[inix][1])) / 2
                xx1 = tempp1[inix][1] - tempp
                xx3 = tempp1[inix + 1][1] + tempp12

            tempp2.append((tempp1[inix][0], tempp1[inix][2], tempp1[inix][3], xx1, tempp1[inix][1], xx3))
        newml[q1] = tempp2

    # xxx = uav_deal_ang(uav_pos[ss],[init_position_matrix[en]])
    mm1 = []
    for a1 in alnodess:
        ff = int(a1)

        # node_time.append(ntt[0][0])
        s = -1
        k1 = uav_deal_ang(locm, [init_position_matrix[ff]])
        cd1 = 0
        if (a1 in node_num):
            # node_inf.append((node[0], theta, dn, k, uav))
            dnq = math.sqrt(pow(init_position_matrix[ff][0, 1] - uav_pos[jj][0], 2) + pow(
                init_position_matrix[ff][0, 2] - uav_pos[jj][1], 2))
            theta1 = math.atan(398.5 / dnq)
            cd1 = 1

        else:
            dnq = math.sqrt(pow(init_position_matrix[ff][0, 1] - uav_pos[jj][0], 2) + pow(
                init_position_matrix[ff][0, 2] - uav_pos[jj][1], 2))
            theta1 = math.atan(398.5 / dnq)
            if (dnq < 990.2):
                cd1 = 0
            else:
                cd1 = -1

        temp1 = []
        for b1 in range(len(zob_map[a1])):
            s = s + 1

            if (zob_map[a1][0][4] > k1[0][1]):
                t1 = (zob_map[a1][0][4] - k1[0][1]) / 2
                x1 = k1[0][1] - t1
                x3 = k1[0][1] + t1
                if (k1[0][1] == 0):
                    x1 = 0
                elif (x1 < 0):
                    x1 = 0
                temp1.append((len(uav_pos)-1, theta1, cd1, x1, k1[0][1], x3))
                temp1.append((zob_map[a1][0][0], zob_map[a1][0][1], zob_map[a1][0][2], zob_map[a1][0][3] - t1,
                              zob_map[a1][0][4], zob_map[a1][0][5]))
                for c1 in range(len(zob_map[a1]) - 1):
                    temp1.append(zob_map[a1][c1 + 1])
                newml[a1] = temp1
                break
            elif (zob_map[a1][len(zob_map[a1]) - 1][4] < k1[0][1]):
                t1 = (k1[0][1] - zob_map[a1][len(zob_map[a1]) - 1][4]) / 2
                x1 = k1[0][1] - t1
                x3 = k1[0][1] + t1
                for c1 in range(len(zob_map[a1]) - 1):
                    temp1.append(zob_map[a1][c1])
                temp1.append((zob_map[a1][len(zob_map[a1]) - 1][0], zob_map[a1][len(zob_map[a1]) - 1][1],
                              zob_map[a1][len(zob_map[a1]) - 1][2], zob_map[a1][len(zob_map[a1]) - 1][3],
                              zob_map[a1][len(zob_map[a1]) - 1][4], zob_map[a1][len(zob_map[a1]) - 1][4] + t1))
                if (k1[0][1] == math.pi * 2 or x3 >= math.pi * 2):
                    x3 = math.pi * 2
                temp1.append((len(uav_pos)-1, theta1, cd1, x1, k1[0][1], x3))
                newml[a1] = temp1
                break
            elif (zob_map[a1][s][4] < k1[0][1] and zob_map[a1][s + 1][4] > k1[0][1]):
                for c1 in range(len(zob_map[a1])):
                    if (c1 == s):
                        t1 = (k1[0][1] - zob_map[a1][s][4]) / 2
                        t2 = (zob_map[a1][s + 1][4] - k1[0][1]) / 2
                        x1 = k1[0][1] - t1
                        x3 = k1[0][1] + t2
                        temp1.append((zob_map[a1][s][0], zob_map[a1][s][1],
                                      zob_map[a1][s][2], zob_map[a1][s][3],
                                      zob_map[a1][s][4], zob_map[a1][s][4] + t1))

                        temp1.append((len(uav_pos)-1, theta1,cd1, x1, k1[0][1], x3))
                        temp1.append((zob_map[a1][s + 1][0], zob_map[a1][s + 1][1],
                                      zob_map[a1][s + 1][2], zob_map[a1][s + 1][4] - t2,
                                      zob_map[a1][s + 1][4], zob_map[a1][s + 1][5]))
                    elif (c1 > s + 1 or c1 < s):
                        temp1.append(zob_map[a1][c1])
                newml[a1] = temp1
                break

    '''
    if (len(shiji) >= 45):
        solveH(locm[0], locm[1], alnodess, init_position_matrix, newml, shiji)
    '''
    print("------------------------split line------------------------")
    demo = alfour.PSO(alnodess, init_position_matrix, newml, ctime)
    demo.fanwei()
    demo.initial()
    demo.solving(15)
    result = demo.returnbest()
    result.sort(key=takeFirst, reverse=True)
    # hij = solveH(result[0][1], result[0][2], alnodess, init_position_matrix, newml)
    locm = [result[0][1], result[0][2], 400]
    print(result[0])
    if (ctime > -6):
        adaptive_dep(uav_pos, newml, init_position_matrix, locm, alnodess, ctime, suoyou)
    else:
        node_inti, node_inf2, node_num, shiji, alloss = bestpp(locm, init_position_matrix)
        print('范围内', len(node_num))
        print('实际符合：', len(shiji))
        print('laaffbf')
        suoyou.append((locm, len(node_num), node_num))
        mdt = []
        for duiz in suoyou:
            mdt.append(duiz[1])
        mdt.sort(reverse=True)
        mdta = []
        for ala in suoyou:
            if(ala[1]>=mdt[0]):
                mdta.append(ala)
        aiya = []
        if(len(mdta)>1):
            for ss in mdta:
                tdis = 0
                for kk in ss[2]:
                    tdis = tdis + pow(ss[0][0]-init_position_matrix[int(kk)][0,1],2) + pow(ss[0][1]-init_position_matrix[int(kk)][0,2],2)
                aiya.append((ss[0],tdis,ss[2]))
            aiya.sort(key=takeSecond)
            H_ite(aiya[0][0],aiya[0][2], alnodess, init_position_matrix, newml, 1)
        else:
            H_ite(mdta[0][0],mdta[0][2], alnodess, init_position_matrix, newml,1)

def H_ite (xyz,minnode, alnodes, init_position_matrix, zob_map,time):
    print('liuyuheng')
    yuanshi = len(minnode)
    mubiao = []
    mubiao.append((xyz,yuanshi,yuanshi))
    if(time==1):
        xxx1 = uav_deal_ang(xyz, init_position_matrix)

        tiao = 0
        haiya = []
        for tatq in alnodes:
            # print('有点问题:',zob_map[int(tatq)])
            for tat1 in zob_map[int(tatq)]:
                if (xxx1[int(tatq)][1] >= tat1[3] and xxx1[int(tatq)][1] <= tat1[5]):
                    if(tat1[2] == -1):
                        tiao = 1
                        tee = pow(pow(1068, 2) - pow(init_position_matrix[int(tatq)][0, 1] - xyz[0], 2) - pow(
                            init_position_matrix[int(tatq)][0, 2] - xyz[1], 2),0.5)
                        haiya.append((tatq,tee))
                        break
        haiya.sort(key=takeSecond)

        jjs = 0
        if(haiya!=[]):
            print('laishishi')
            while jjs==0:
                xyz = (xyz[0], xyz[1], xyz[2] - 10)
                node_inti, node_inf2, node_num, shiji, alloss = bestpp(xyz, init_position_matrix)
                if(len(shiji)>=yuanshi):
                    mubiao.append((xyz,len(shiji),shiji))
                if(xyz[2]<haiya[0][1]):
                    jjs = 1
                elif(len(shiji)+len(haiya)<=yuanshi):
                    jjs = 1
            mubiao.sort(key=takeSecond,reverse=True)
            xyz = mubiao[0][0]
        else:
            print('eiieie')
            while jjs==0:
                xyz = (xyz[0], xyz[1], xyz[2] - 10)
                node_inti, node_inf2, node_num, shiji, alloss = bestpp(xyz, init_position_matrix)
                if(len(shiji)<yuanshi):
                    jjs = 1
                    xyz = (xyz[0], xyz[1], xyz[2] + 10)
                if(xyz[2]<300):
                    jjs = 1
                    xyz = (xyz[0], xyz[1], xyz[2] + 10)

        print('最终的位置', xyz)
        node_inti, node_inf2, node_num, shiji, alloss = bestpp(xyz, init_position_matrix)
        print(len(shiji))
        #engc(xyz, init_position_matrix)



def engc(xyz,init_position_matrix):
    print('finial')
    node_inti, node_inf2, node_num, shiji, alloss = bestpp(xyz, init_position_matrix)
    ptlist = []
    x = Symbol('x')
    for jisuan in node_inf2:
        #  node_inf.append((node[0], theta, dn, 1, uav))

        PL = 20 * log(jisuan[2]) + 20 * math.log10(2) + 92.4
        pt = 10 ** ((-70 + x + PL) / 10)
        ft = exp(-x ** 2 / 32) / ((2 * pi) ** 0.5 * 4)
        ca = 1-integrate(ft, (x, float('-inf'), -x))
        Nj = pt / (ca)
        m = diff(Nj, x)
        ktt = 0
        s = -1
        while ktt==0:
            s = s+0.005
            y_1 = m.evalf(subs={x: s})
            if(y_1>=0):
                ptlist.append((-70 + s + PL,jisuan[0]))
                print(s)
                ktt = -1



def takeSecond(elem):
    return elem[1]


def takeFirst(elem):
    return elem[0]


time1 = time.time()
node_list = []
com_node_list = []
bian, tnum, ditu, A = pre_grid(1)

# movement_matrix, init_position_matrix = Gm.get_position('ggnodes121.tcl')
# movement_matrix, init_position_matrix = Gm.get_position('nodes-type2.tcl')
# movement_matrix, init_position_matrix = Gm.get_position('ggnodes11.tcl')
# movement_matrix, init_position_matrix = Gm.get_position('ggnodes.tcl')

#movement_matrix, init_position_matrix = Gm.get_position('2020-8-16-52-2.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-3-10-type2.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-3-11-type2.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-7-24-52type2.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-13-46-node2.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-13-51-node2.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-16-52-2.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-7-26-t.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-23-50-cc.tcl')

#movement_matrix, init_position_matrix = Gm.get_position('2020-8-25-x-47.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-25-x-49.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-25-x-57.tcl')

#movement_matrix, init_position_matrix = Gm.get_position('2020-8-27-d-55.tcl')
movement_matrix, init_position_matrix = Gm.get_position('2020-8-25-x-t-50.tcl')
node_num = init_position_matrix.shape[0]
print('num:', node_num)

nn = [500, 500, 400]
xxx = []
temp = []
all = []
enenti = 0
node_inti, node_inf, node_num, shiji = broadcast1(nn, init_position_matrix, 0, -70, 2, 1000)
print('范围内', len(node_num))
print('实际符合：', len(shiji))
k = 0
for dtemp in A:
    if (dtemp[0] != 500 and dtemp[1] != 500):
        tedis = math.sqrt(pow(dtemp[0] - nn[0], 2) + pow(dtemp[1] - nn[1], 2)) / 1000
        if (tedis < 0.5):
            dddd = kenengxing2(tedis)
            temp4 = dtemp[2] * dddd
            A[k] = [dtemp[0], dtemp[1], temp4]
            if (A[k][2] < 0.02):
                A[k] = [500, 500]
    k = k + 1

temp.append((nn,len(node_num),copy.deepcopy(node_num)))
xxx.append(node_inf)
print(temp)

neib = []
UAV_Information_collection(nn, init_position_matrix, A, neib, enenti, temp, all, xxx)

time2 = time.time()
print('END')
print(time2 - time1)

