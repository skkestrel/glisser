# -*- coding: UTF-8 -*-
import numpy as np
from random import *
from elements import *


def generateRandomStates(arange = [50,100], num = 1000, eles = [], states = []):
    for i in range(num):
        a = np.random.uniform(arange[0], arange[1])
        e = np.random.uniform(0,0.4)
        # q = np.random.uniform(28, 34)
        # e = 1 - q/a
        I = np.random.uniform(0, 20)
        # a = a0+step*i
        ome = np.random.uniform(0, 360)
        Ome = np.random.uniform(0, 360)
        f = np.random.uniform(0, 360)
        # f = 180
        ele = [a, e, I, ome, Ome, f]
        print(ele)
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


def writeGlisseSMALLIN(states, starting = 1, filename="GLISSER-small.in"):
    file = open(filename, "w+")
    if starting == 1:
        file.write(str(len(states))+"\n")
    for i, state in zip(np.arange(len(states)), states):
        for dat in state[:3]:
            file.write(" {: 2.15f}".format(dat))
        file.write("\n")
        for dat in state[3:]:
            file.write(" {: 2.15f}".format(dat))
        file.write("\n")
        file.write(" "+str(i+starting)+" 0 0\n")
    file.close()

codename = "Resonant-Test"

eles, states = generateRandomStates([20,100], 1000)
# eles, states = generateRandomStates([47.5,47.9], 1000, eles, states)
# eles, states = generateRandomStates([55.2,55.6], 1000, eles, states)
# eles, states = generateRandomStates([62.3,62.7], 1000, eles, states)
# writeSwiftTPIN(states)
# writeMercurySMALLIN(eles)
writeGlisseSMALLIN(states, filename="GLISSER-small-REBOUND.in".format(codename))


# ele = [ -20.749197946439399, -19.795374429795853, -39.180417442715310,
#   0.001683063894229  ,0.001098110929503 ,-0.001423157829402]
# print(ele)
# state = c2e(ele, mu_glisse)
# print(state)
