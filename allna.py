# coding: utf-8
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import hull

import math
import random
import matplotlib.pyplot as plt


#500米高
c1 = 2  # 学习因子                          # Learning coefficient
c2 = 2


def fitness(x1, x2, x3):  # 适应度函数                            # Fitness Function
    return math.floor((2 * x1 ** 2 - 3 * x2 ** 2 - 4 * x1 + 5 * x2 + x3) * 100) / 100


class PSO:
    def __init__(self, alnodes, init_position_matrix, map, tim):
        self.pop_size = 100 # 粒子群个体数                         # the number of the instance in Particle swarm
        self.dim = 2  # 变量数                               # the number of variables
        self.omega = 0.6  # 惯性因子                             # Inertia factor
        self.x_max = 900
        self.x_min = 100
        self.y_min = 0
        self.y_max = 0
        self.vx_max = (self.x_max - self.x_min) * 0.02
        self.vx_min = -(self.x_max - self.x_min) * 0.02
        self.vy_max = (self.y_max - self.y_min) * 0.02
        self.vy_min = -(self.y_max - self.y_min) * 0.02
        self.position = [[]]  # 记录当前粒子位置                      # record the current position of each particle
        self.speed = [[]]  # 记录当前粒子速度                      # record the moving speed
        self.best_value = [[]]  # 记录全局最优                          # record the global optimal
        self.value = [[]]  # 记录当前值                            # record the current fitness and position
        self.alnodes = alnodes
        self.init_position_matrix = init_position_matrix
        self.map = map
        self.tim = tim

    def initial(self):
        for i in range(self.pop_size):
            "第一轮初始化"
            "first round of initialization"
            x = []
            v = []
            x.append(math.floor(random.uniform(self.x_min, self.x_max) * 100) / 100)
            v.append(math.floor(random.uniform(self.vx_min, self.vx_max) * 100) / 100)
            x.append(math.floor(random.uniform(self.y_min, self.y_max) * 100) / 100)
            v.append(math.floor(random.uniform(self.vy_min, self.vy_max) * 100) / 100)
            self.position.append(x)
            self.speed.append(v)
            # self.value.append((fitness(x[0],x[1],x[2]))
            self.value.append(((self.objfun([x[0], x[1]])), x[0], x[1]))
        self.value = self.value[1:]
        # -------------------选择取最大值/最小值-------------------#
        "choose to get the max or the min value"
        index = self.value.index(max(self.value))
        # index = self.value.index(min(self.value))
        # -------------------选择取最大值/最小值-------------------#
        self.best_value.append((self.value[index][0], self.position[index][0], self.position[index][1]))
        self.best_value = self.best_value[1:]
        self.position = self.position[1:]
        self.speed = self.speed[1:]
        # print("the population and fitness after initialization:")
        # print("position :",self.position)
        # print("value:",self.value)
        # print("best value:",self.best_value)

    def fanwei(self):
        mm1 = []
        for a1 in self.alnodes:
            ff = int(a1)
            mm1.append((self.init_position_matrix[ff][0, 1], self.init_position_matrix[ff][0, 2]))
        mm1.append((self.init_position_matrix[int(self.alnodes[0])][0, 1],
                    self.init_position_matrix[int(self.alnodes[0])][0, 2]))
        oo = []
        fanwei = hull.convex(mm1)
        list.sort(fanwei, key=self.takeFirst)
        self.x_min = fanwei[0][0]
        self.x_max = fanwei[len(fanwei) - 1][0]
        list.sort(fanwei, key=self.takeSecond)
        self.y_min = fanwei[0][1]
        self.y_max = fanwei[len(fanwei) - 1][1]

    def takeSecond(self, elem):
        return elem[1]

    def takeFirst(self, elem):
        return elem[0]

    def solving(self, times):
        for i in range(times):
            # print(self.value)
            # -------------------选择取最大值/最小值-------------------#
            "choose to get the max or the min value"
            pbest = self.value[self.value.index(max(self.value))]
            gbest = self.best_value[self.best_value.index(max(self.best_value))]
            # pbest = self.value[self.value.index(min(self.value))]
            # gbest = self.best_value[self.best_value.index(min(self.best_value))]
            # -------------------选择取最大值/最小值-------------------#
            # print("pbest:",pbest)
            # print("gbest",gbest)
            for j in range(self.pop_size):
                x = []
                v = []
                for k in range(self.dim):
                    v.append(math.floor((self.omega * self.speed[j][k] + c1 * random.uniform(0, 1) * (
                                pbest[1 + k] - self.position[j][k]) + c2 * random.uniform(0, 1) * (
                                                     gbest[1 + k] - self.position[j][k])) * 100) / 100)
                    x.append(math.floor((self.position[j][k] + self.speed[j][k]) * 100) / 100)
                    "将位置和速度限制在规定范围内"
                    "restrict the position and the speed"
                    if (k == 0):
                        if (v[k] < self.vx_min):
                            v[k] = self.vx_min
                        if (v[k] > self.vx_max):
                            v[k] = self.vx_max
                        if (x[k] < self.x_min):
                            x[k] = self.x_min
                        if (x[k] > self.x_max):
                            x[k] = self.x_max
                    else:
                        if (v[k] < self.vy_min):
                            v[k] = self.vy_min
                        if (v[k] > self.vy_max):
                            v[k] = self.vy_max
                        if (x[k] < self.y_min):
                            x[k] = self.y_min
                        if (x[k] > self.y_max):
                            x[k] = self.y_max
                "数据更新"
                "updating"
                self.position[j] = x
                self.speed[j] = v
                self.value[j] = (
                self.objfun((self.position[j][0], self.position[j][1])), self.position[j][0], self.position[j][1])
            # -------------------选择取最大值/最小值-------------------#
            "choose to get the max or the min value"
            index = self.value.index(max(self.value))
            # index = self.value.index(min(self.value))
            # -------------------选择取最大值/最小值-------------------#
            "获得全局最优"
            "get the global optimal"
            self.best_value.append((self.value[index][0], self.position[index][0], self.position[index][1]))

    def returnbest(self):
        return self.best_value

    def uav_deal_ang(self, uav, inf):

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

    def objfun(self, post):
        # 571
        sumcl = 0
        ll = 0
        l = 0

        cc1 = 10 * (10 - self.tim)
        for max in self.map:

            if (max in self.alnodes):
                jii = 0
                l = l + 1
                temij = self.uav_deal_ang(post, self.init_position_matrix[max])
                dis = math.sqrt(pow(post[0] - self.init_position_matrix[max][0, 1], 2) + pow(
                    post[1] - self.init_position_matrix[max][0, 2], 2))
                tat = math.degrees(math.atan(698.5 / dis))
                plos1 = (120 - (120 / (1 + pow(tat / 24.3, 1.229)))) / 100

                for suan in range(len(self.map[max])-1):
                    if (temij[0][1] > self.map[max][suan][3] and temij[0][1] <= self.map[max][suan][5]and dis< 806.6):
                        if (self.map[max][suan][2] == 1):
                            sumcl = sumcl + 1
                            break

                        elif (self.map[max][suan][2] == 0 and tat > self.map[max][suan][1]):
                            ll = ll +  plos1
                            jii = 1
                            break


                        elif (self.map[max][suan][2] == -1 and tat > self.map[max][suan][1]):
                            ll = ll +  plos1
                            jii = 1
                            break





        if (self.tim >= 6):
            return sumcl + ll

        else:
            return sumcl


        #return sumcl


'''
if __name__ == '__main__':
    demo = PSO()
    demo.initial()
    demo.solving(100)
    result = demo.returnbest()
    for i in result:
        print("value:",i[0],"x1:",i[1],"x2:",i[2])
        
    if (tat >= self.map[max][suan][1]):
        sumcl = sumcl + 1
    else:
        sumcl = sumcl + 1 * plos1

'''