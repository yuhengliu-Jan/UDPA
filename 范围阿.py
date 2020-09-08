import time
import Get_Move as Gm
from pyscipopt import Model, quicksum, multidict, SCIP_PARAMSETTING, exp, log, sqrt


#movement_matrix, init_position_matrix = Gm.get_position('2020-8-25-x-47.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-25-x-49.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-25-x-57.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-27-d-55.tcl')
movement_matrix, init_position_matrix = Gm.get_position('2020-8-25-x-t-50.tcl')

#movement_matrix, init_position_matrix = Gm.get_position('2020-8-13-46-node2.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-13-51-node2.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-16-52-2.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-7-26-t.tcl')
#movement_matrix, init_position_matrix = Gm.get_position('2020-8-23-50-cc.tcl')

print(len(init_position_matrix))

model = Model()  # model name is optional
cl, cn, ui = {}, {}, {}
themx, thetav, theta, hij = {}, {}, {}, {}
for cj in range(len(init_position_matrix)):
    ui[cj] = model.addVar(vtype="B", name="bij(%s)" % (cj))


#hij = model.addVar(vtype="C", name="hij(%s)" )
xij = model.addVar(vtype="C", name="xij(%s)")
yij = model.addVar(vtype="C", name="yij(%s)")
'''
300 107637.56687421931
400 191355.67444305654
500 298993.2413172759
600 430550.26749687723
700 586026.7529818608
mijichengs
300 45386.61010972501
400 80687.30686173333
500 126073.91697145837
600 181546.44043890003
700 247104.8772640584
mijicheng hai gao
300 6002.152755492015
400 10670.49378754136
500 16672.64654303337
600 24008.61102196806
700 32678.38722434541
'''
for kk in range(len(init_position_matrix)):
    model.addCons((init_position_matrix[kk][0,1]-xij)**2 + (init_position_matrix[kk][0,2]-yij)**2 <= 107637.56687421931 + 1000000000*(1-ui[kk]) )

model.setObjective(sum(ui[i] for i in ui),sense="maximize")
#model.setObjective(sum(ui[i]*(init_position_matrix[i][0,1]-xij)**2 + (init_position_matrix[i][0,2]-yij)**2 for i in ui),sense="minimize")
model.optimize()

objSet = bool(model.getObjective().terms.keys())
print("* Is objective set? %s" % objSet)
#print(',,,',model.getVal(xij))
#print(',,,',model.getVal(yij))
for kk in model.getVars():
    print("%s: %d" % (kk, round(model.getVal(kk))))

print(';;;;')

