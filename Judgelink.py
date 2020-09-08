#Python3.6
import numpy as np
from numpy.linalg import solve
import math
import re
import matplotlib.pyplot as plt

class point(): #定义类
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z

def cross(p1,p2,p3):#跨立实验
    x1=p2.x-p1.x
    y1=p2.y-p1.y
    x2=p3.x-p1.x
    y2=p3.y-p1.y
    return x1*y2-x2*y1

def IsIntersec(p1,p2,p3,p4): #判断两线段是否相交

    #快速排斥，以l1、l2为对角线的矩形必相交，否则两线段不相交
    if(max(p1.x,p2.x)>=min(p3.x,p4.x)    #矩形1最右端大于矩形2最左端
    and max(p3.x,p4.x)>=min(p1.x,p2.x)   #矩形2最右端大于矩形最左端
    and max(p1.y,p2.y)>=min(p3.y,p4.y)   #矩形1最高端大于矩形最低端
    and max(p3.y,p4.y)>=min(p1.y,p2.y)): #矩形2最高端大于矩形最低端

    #若通过快速排斥则进行跨立实验
        if(cross(p1,p2,p3)*cross(p1,p2,p4)<=0
           and cross(p3,p4,p1)*cross(p3,p4,p2)<=0):
            D=1   #相交了
        else:
            D=0
    else:
        D=0
    return D


def Lin_eq(x1, x2):  #通过两点计算二维直线方程得到y-kx=b的形式
  if (x2.x - x1.x==0):
      return [0,1,x1.x]
  elif(x2.y - x1.y==0):
      return [1, 0, x1.y]
  else:
      k = (x2.y - x1.y) / (x2.x - x1.x)
      remnum = x1.y - k*x1.x
      k = 0 - k
      return [1, k, remnum]

def Solu_eq(z1,z2,z3,z4):   #如果二维两直线相交计算交点
    node1 = Lin_eq(z1,z2)
    node2 = Lin_eq(z3,z4)
    A1 = np.mat([[node1[1], node1[0]], [node2[1], node2[0]]])
    B1 = np.mat([node1[2], node2[2]]).T
    xs = solve(A1, B1)
    i = point(float(xs[0][0]),float(xs[1][0]),0)
    return i

def Seek_alt(n1,n2,n3): #判断三维空间中两点之间是否有建筑物
    t1 = (n2.z - n1.z) / (n2.x - n1.x)
    h =  (n3.x - n2.x)*t1
    if (h>n3.z): D = 0;
    else: D = 1 #有建筑物
    return D

def Judge_h(u1,u2,g1,g2):    #两对节点
    Ju1 = IsIntersec(u1,u2,g1,g2) # 1相交了
    if(Ju1 == 1):
       Inter = Solu_eq(u1,u2,g1,g2)
       n3 = point(0,0,0)
       n3.x = Inter.x
       n3.y = Inter.y
       n3.z = g1.z
       pa = Seek_alt(u1,u2,n3)
       return pa
    else:
        return 0


def J_broadcast(a1, a2):
    # print('我还没写')
    bu = 1 #1没相交0相交了
    if (a1.x==a2.x and a1.y == a2.y):
        return bu
    build_list = get_build('map_inf.tcl')
    for build in build_list:
        A1 = point(build[0],build[1],build[8])
        A2 = point(build[2], build[3], build[8])
        A3 = point(build[4], build[5], build[8])
        A4 = point(build[6], build[7], build[8])
        temp1 = Judge_h(a1,a2,A1,A2)
        temp2 = Judge_h(a1, a2, A2, A3)
        temp3 = Judge_h(a1, a2, A3, A4)
        temp4 = Judge_h(a1, a2, A1, A4)
        if (temp1==1 or temp2==1 or temp3==1 or temp4==1):
            # print(build[8])
            bu = 0
            break


    return bu


def gns_dis(position_list):
    # 地面节点之间的距离列表
    list = []
    sum = []
    k = 0
    node_num = position_list.shape[0]
    # print(node_num)
    for node in position_list:
        aimn = point(node[0,1],node[0,2],node[0,3])
        list.append(k)

        for i in range(node_num):
            if (i == k):
                dis = 0
                list.append(dis)
            else:
                #for j in position_list:
                temp = pow((aimn.x-position_list[i][0,1]), 2 ) +  pow((aimn.y-position_list[i][0,2]), 2 ) +  pow((aimn.z-position_list[i][0,3]), 2 )
                dis = math.sqrt(temp)
                # print(dis)
                list.append(dis)

        #list.clear()
        #print(list)
        c = list.copy()
        sum.append(c)
        k = k + 1
        #print(sum)
        list.clear()

    return sum


def Ag_dis(uav_position, gn_position):
    temp = []
    temp1 = []
    n_num = gn_position.shape[0]


def get_build(mobile_file_path):
    build_list = []
    with open(mobile_file_path, 'r') as f:
        for line in f :
            line_list = re.split('[\s]', line)
            # build_list.append((int(float(line_list[0])), int(float(line_list[1]))))
            build_list.append((float(line_list[3]), float(line_list[5]), float(line_list[3]), float(line_list[6]), float(line_list[4]), float(line_list[5]), float(line_list[4]), float(line_list[6]), float(line_list[7])))
            # print(line_list[7])
        # print(build_list)
    return build_list



'''
node1 = point(500 , 500 , 100)
node3 = point(988.8876609269237,124.12661542986869,100)

node2 = point(988.8876609269237,124.12661542986869,1.5)


hh = J_broadcast(node1,node2)
print(hh)

m1 = point(0,0,0)
m2 = point(0,1,0)
m3 = point(1,0,0)
m4 = point(1,1,0)


Ju12 = IsIntersec(m1,m2,m3,m4)
print(Ju12)'''



'''
mm = point(50,50,3)
ma = point(50,100,100)
ss = J_broadcast(ma,mm)
print(ss)'''