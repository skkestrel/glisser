# -*- coding: UTF-8 -*-
from re import M
import numpy as np
import rebound
import orbitplotkit as opk
import invariant_plane as ivp
import pandas as pd
import struct
import time
import sys
import random
import os

print(rebound.__version__)

co = 0.00029592338593516655
co2 = 365.25/(2*np.pi)
date = "2022-07-01 00:00"
output_folder = os.path.join(os.getcwd())

print(output_folder)


def genSim():

    sim = rebound.Simulation()
    sim.add("Sun", date=date)
    sim.add("Mercury", date=date)
    sim.add("Venus", date=date)
    sim.add("Earth", date=date)
    sim.add("Mars", date=date)
    sim.add("Jupiter", date=date)
    sim.add("Saturn", date=date)
    sim.add("Uranus", date=date)
    sim.add("Neptune", date=date)

    sim.save("sim_eight_planets_20220701.bin")
    NUM = sim.N
    jpl_data = pd.read_csv("sbdb_query_results_8.23.csv")
    jpl_data = jpl_data[jpl_data['condition_code'] <= 5]
    # jpl_data = jpl_data[jpl_data['q'] > 36]
    jpl = jpl_data[jpl_data['a'] > 50].copy()

    # print(jpl)

    for i in np.arange(jpl[jpl.columns[0]].count()):
        sim.add(jpl.iloc[i]['full_name'].split('(', 1)[1].split(')')[0], m=0)
    # sim.add('2002 CP154', m=0)
    orbits = sim.calculate_orbits(primary=sim.calculate_com())
    a_bary, e_bary, inc_bary, q_bary = [], [], [], []
    for orbit in orbits[NUM-1:]:
        print(orbit.a)
        a_bary.append(orbit.a)
        e_bary.append(orbit.e)
        inc_bary.append(orbit.inc)
        q_bary.append(orbit.a*(1-orbit.e))
    jpl['a_bary'] = a_bary
    jpl['e_bary'] = e_bary
    jpl['inc_bary'] = inc_bary
    jpl['q_bary'] = q_bary
    jpl.to_csv("distant_tno_bary_8.23.csv")

    return sim


sim = genSim()

# genBIG(sim, os.path.join(output_folder, "big-JSUNT-{0}".format(M)))

# genTestHist(sim, M, 200, 40180000000, 1000, 100, True)

# genHist("plhist_q=34_M=40", 200, 73200000000, 1000, True)
