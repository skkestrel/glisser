# -*- coding: UTF-8 -*-
import numpy as np
from random import *
from elements import *


def generateRandomStates():
    eles = []
    states = []
    for i in range(num):
        a = np.random.uniform(4, 50)
        # q = 50
        e = np.random.uniform(0, 0.01)
        I = np.random.uniform(80, 100)
        # a = a0+step*i
        ome = np.random.uniform(0, 360)
        Ome = np.random.uniform(0, 360)
        f = np.random.uniform(0, 360)
        # f = 180
        ele = [a, e, I, ome, Ome, f]
        eles.append(ele)
        state = e2c(ele, mu_glisse)
        states.append(state)
    return eles, states


def writeSwiftTPIN(states):
    file = open("tp_swift.in", "w+")
    file.write(str(num)+"\n")
    for state in states:
        for dat in state[:3]:
            file.write(" {: 2.15f}".format(dat))
        file.write("\n")
        for dat in state[3:]:
            file.write(" {: 2.15f}".format(dat))
        file.write("\n")
        file.write("  0 0 0 0 0 0 0 0 0 0 0 0 0\n")
        file.write("  0.0d0 0.0d0 0.0d0 0.0d0 0.0d0\n")
        file.write("  0.0d0 0.0d0 0.0d0 0.0d0 0.0d0\n")
        file.write("  0.0d0 0.0d0 0.0d0\n")
    file.close()


def writeMercurySMALLIN(eles):
    file = open("yh_small.in", "w+")
    file.write(
        ")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\n")
    file.write(
        ") Lines beginning with `)' are ignored.\n")
    file.write(
        ")-------------------------------------------------------------------\n")
    file.write(
        " style(Cartesian, Asteroidal, Cometary)=Asteroidal\n")
    file.write(
        ")------------------------d or D is not matter--0d0 is possible too--\n")
    for i, ele in zip(np.arange(num), eles):
        ele[5] = np.rad2deg(true2mean(np.deg2rad(ele[5]), ele[1]))
        file.write(" {0:d}    ep={1:.1f}d0\n".format(i+101, 0))
        for dat in ele[:3]:
            file.write(" {: 2.15f}".format(dat))
        file.write("\n")
        for dat in ele[3:]:
            file.write(" {: 2.15f}".format(dat))
        file.write("\n")
        file.write(" 0.d0 0.d0 0.d0\n")
    file.close()


def writeGlisseSMALLIN(states):
    file = open("GLISSER-small-{0}.in".format(codename), "w+")
    file.write(str(len(states))+"\n")
    for i, state in zip(np.arange(num), states):
        for dat in state[:3]:
            file.write(" {: 2.15f}".format(dat))
        file.write("\n")
        for dat in state[3:]:
            file.write(" {: 2.15f}".format(dat))
        file.write("\n")
        file.write(" "+str(i+1)+" 0 0\n")
    file.close()

codename = "Polar"
num = 80000

# eles, states = generateRandomStates()
# writeSwiftTPIN(states)
# writeMercurySMALLIN(eles)
# writeGlisseSMALLIN(states)


ele = [ -20.749197946439399, -19.795374429795853, -39.180417442715310,
  0.001683063894229  ,0.001098110929503 ,-0.001423157829402]
print(ele)
state = c2e(ele, mu_glisse)
print(state)
