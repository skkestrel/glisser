#!/usr/bin/env python3
# import numpy as np
# import matplotlib.pyplot as plt
import sys
import time
import os

CONFIG_OLD = """
Initial-Time 0
Time-Step 100
Final-Time 100000000
Time-Block-Size 100
Resolve-Encounters 1
Planet-History-File benchmark/outer-planets-out/plhist.out
Read-Planet-History 1
Swift-Path swift/main/swift_glisse_ext
Track-Interval 1
Encounter-RH-Factor 3.5
Cull-Radius 1
Log-Interval 100
Dump-Interval 20000
Input-File benchmark/outer-particles/state.in
Swift-Process-Count 8
Swift-Process-Min-Particle-Count 10
"""

CONFIG = """
Initial-Time 0
Time-Step 100
Final-Time 20000000
Time-Block-Size 100
Resolve-Encounters 1
Planet-History-File benchmark/outer-planets-out/plhist.out
Read-Planet-History 1
Swift-Path swift/main/swift_glisse_ext
Track-Interval 0
Resync-Interval 1
Encounter-RH-Factor 1
Cull-Radius 1
Log-Interval 500
Dump-Interval 0
Input-File benchmark/circular-particles/state.in
Swift-Process-Count 8
Swift-Process-Min-Particle-Count 10
"""

os.system('rm benchmark/timelog.out')
os.system('rm -r benchmark/circular-particles-out-*')

num = 100
numlist = []
while(num <= 140000):
    numlist.append(num)
    if(num < 1000):
        num += 100
    elif(num < 50000):
        num += 1000
    elif(num <= 200000):
        num += 2000
print(numlist)

for lps in numlist:
    ENDLINE1 = "Limit-Particle-Count " + str(lps)
    ENDLINE2 = "Output-Folder benchmark/circular-particles-out-"+ str(lps)
    CONFIG_NEW = CONFIG + ENDLINE1 + '\n' + ENDLINE2
    print(CONFIG_NEW)

    os.system('rm benchmark/config.in')
    fin = open("benchmark/config.in", "w")
    fin.write(CONFIG_NEW)
    fin.close()

    now = time.time()
    os.system('./bin/glisser benchmark/config.in')
    difference = time.time() - now
    print(difference)

    fout = open("benchmark/timelog.out", "a")
    fout.write(str(lps) + " " + str(difference) + "\n")
    fout.close()