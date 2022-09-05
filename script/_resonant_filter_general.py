import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import os
import matplotlib
import orbitplotkit as opk


font = {'weight': 'bold',
        'size': 16}
matplotlib.rc('font', **font)


folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Giant-Synthetic/"
# folder = os.path.join(os.getcwd(), "fast" , "rogue_output", "out-JSUNR-Resonant-E")

output_folder = os.path.join(folder, "pics", "hist")
hist_folder = os.path.join(folder, "reoutput", "hist")
histRESO_folder = os.path.join(folder, "reoutput", "hist_RESO")
snapshot_folder = os.path.join(folder, "reoutput", "snapshots")
snapshotRESO_folder = os.path.join(folder, "reoutput", "snapshots_RESO")
enc_folder = os.path.join(folder, "reoutput", "enc")
encounter_file = "encounter.out"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

if not os.path.exists(enc_folder):
    os.makedirs(enc_folder)

if not os.path.exists(snapshotRESO_folder):
    os.makedirs(snapshotRESO_folder)

if not os.path.exists(histRESO_folder):
    os.makedirs(histRESO_folder)


class Filter:

    def __init__(self, npl, target_pl, npa, cadence, windowA, windowB, shift, filename, isOutputHist):
        self.maxnpl = npl
        self.maxnpa = npa
        self.npl = npl
        self.npa = npa
        self.target_pl = target_pl
        self.cadence = cadence
        self.windowA = windowA
        self.windowB = windowB
        self.shift = shift
        self.length = 0
        self.timestamp = 0
        self.countB = 0
        self.t1, self.t2, self.t3 = self.timestamp, self.timestamp + \
            windowB, self.timestamp+windowA
        self.alist = []
        self.alonglist = []
        self.ameanlist = []
        self.asumlist = []
        self.amaxlist = []
        self.aminlist = []
        self.code = []
        self.RESOlist = []
        self.filename = filename
        self.isOutputHist = isOutputHist
        for i in np.arange(self.maxnpa+1):
            self.alist.append([])
            self.alonglist.append([])
            self.ameanlist.append([])
            self.asumlist.append([])
            self.amaxlist.append([])
            self.aminlist.append([])
            self.code.append([])
            self.RESOlist.append([])
        self.deathlist = []

    def printAlist(self):
        for i, a in enumerate(self.ameanlist):
            print(i, a)

        for i, a in enumerate(self.code):
            print(i, len(a))

    def readFirstSnapshots(self):
        flag = True
        while flag and (self.timestamp <= self.t2):
            print("Reading timestamp: {0}".format(self.timestamp))
            try:
                df_pl = pd.read_csv(os.path.join(snapshot_folder, "planets_{0:d}.txt".format(self.timestamp)),
                                    names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')
                df_pa = pd.read_csv(os.path.join(snapshot_folder, "particles_{0:d}.txt".format(self.timestamp)),
                                    names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')
                for i in np.arange(0, self.maxnpa+1):
                    # for i in np.arange(1, 10):
                    try:
                        if i == 0:
                            a_temp = df_pl.loc[self.target_pl]['a']
                        else:
                            a_temp = df_pa.loc[i]['a']
                        self.alist[i].append(a_temp)
                    except KeyError:
                        self.deathlist.append(i)
                        pass
                self.timestamp += self.cadence
                flag = True
            except FileNotFoundError:
                flag = False

        for i in np.arange(0, self.maxnpa+1):
            if i not in self.deathlist:
                self.asumlist[i] = np.sum(self.alist[i])
                self.amaxlist[i] = np.max(self.alist[i])
                self.aminlist[i] = np.min(self.alist[i])
                self.countB = len(self.alist[i])
                self.ameanlist[i] = self.asumlist[i]/self.countB
                self.alonglist[i] = self.alist[i].copy()

    def readNextSnapshot(self, isFirst):
        reslim = 0.08
        print("Reading timestamp: {0}".format(self.timestamp))
        try:
            df_pl = pd.read_csv(os.path.join(snapshot_folder, "planets_{0:d}.txt".format(self.timestamp)),
                                names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')
            df_pa = pd.read_csv(os.path.join(snapshot_folder, "particles_{0:d}.txt".format(self.timestamp)),
                                names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')
            for i in np.arange(0, self.maxnpa+1):
                if i not in self.deathlist:
                    try:
                        if i == 0:
                            a_temp = df_pl.loc[self.target_pl]['a']
                        else:
                            a_temp = df_pa.loc[i]['a']

                        self.alist[i].append(a_temp)
                        self.alonglist[i].append(a_temp)
                        asum = self.asumlist[i]
                        if not isFirst:
                            self.alonglist[i].pop(0)
                        a0 = self.alist[i].pop(0)
                        asum += (a_temp - a0)

                        self.asumlist[i] = asum
                        self.ameanlist[i] = asum/self.countB

                        self.amaxlist[i] = np.max(self.alonglist[i])
                        self.aminlist[i] = np.min(self.alonglist[i])

                        klist, jlist, a_rlist = opk.genOuterRes(
                            self.ameanlist[0], 49, 600, high1=5, high2=90, order_lim=90)
                        if i != 0:
                            idx = np.argmin(np.abs(self.ameanlist[i]-a_rlist))
                            k, j, a_res = klist[idx], jlist[idx], a_rlist[idx]
                            if not isFirst:
                                self.code[i].pop(0)
                            if abs(self.ameanlist[i] - a_res) < reslim:
                                self.code[i].append(k*1000+j)
                            else:
                                self.code[i].append(0)

                    except KeyError:
                        self.deathlist.append(i)
                        pass
            self.timestamp += self.cadence
            return True
        except FileNotFoundError:
            return False

    def loopToWindowA(self):
        flag = True
        while flag and (self.timestamp <= self.t3):
            self.length = self.timestamp
            flag = self.readNextSnapshot(True)

        output_t = self.timestamp - self.length - \
            self.cadence + self.shift*self.cadence
        for temp_t in np.arange(0, output_t, self.cadence):
            print("Output Gray Snapshot: {0}".format(temp_t))
            df_pa = pd.read_csv(os.path.join(snapshot_folder, "particles_{0:d}.txt".format(temp_t)),
                                names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ', nrows=self.maxnpa).set_index('id')
            df_pa['RESO'] = -2
            df_pa.to_csv(os.path.join(snapshotRESO_folder, self.filename +
                                      "_{0:d}.csv".format(temp_t)))

        df_pa = pd.read_csv(os.path.join(snapshot_folder, "particles_{0:d}.txt".format(output_t)),
                            names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ', nrows=self.maxnpa).set_index('id')
        self.outputSnapshotTable(df_pa, output_t)
        if self.isOutputHist:
            self.outputHistTable(df_pa, output_t, np.arange(1, 101), True)

    def loopToEnd(self):
        flag = True
        while flag:
            flag = self.readNextSnapshot(False)
            output_t = self.timestamp - self.length - \
                self.cadence + self.shift*self.cadence
            df_pa = pd.read_csv(os.path.join(snapshot_folder, "particles_{0:d}.txt".format(output_t)),
                                names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ', nrows=self.maxnpa).set_index('id')
            self.outputSnapshotTable(df_pa, output_t)

            if self.isOutputHist:
                self.outputHistTable(df_pa, output_t, np.arange(1, 101), False)

    def outputSnapshotTable(self, df_pa, output_t):
        scalim = 0.0375
        print("Output Snapshot: {0}".format(output_t))
        RESO_list = []
        idlist = df_pa.index.values.tolist()
        for i in idlist:
            # print(self.code[i])
            count = len(self.code[i])
            counts = np.bincount(self.code[i])
            code = np.argmax(counts)
            num = self.code[i].count(code)

            if self.amaxlist[i] - self.aminlist[i] > self.ameanlist[i] * scalim:
                RESO = -1
            elif code == -1:
                RESO = -1
            elif code == 0:
                RESO = 0
            elif code > 0:
                if num/count > 0.5:
                    RESO = code
                else:
                    RESO = 0
            RESO_list.append(RESO)
            self.RESOlist[i] = RESO
            print("{0} {1:.2f} {2} {3:.2f} {4}".format(i, num/count,
                  code, self.amaxlist[i] - self.aminlist[i], RESO))

        df_pa['RESO'] = RESO_list
        df_pa.to_csv(os.path.join(snapshotRESO_folder, self.filename +
                     "_{0:d}.csv".format(output_t)))

    def outputHistTable(self, df_pa, output_t, pa_list, isFirst):
        print("Output Hist: {0}".format(output_t))
        for i in pa_list:
            if i not in self.deathlist:
                pa = df_pa.loc[i]
                if isFirst:
                    output = open(histRESO_folder +
                                  "/pa_{0:d}.csv".format(i), "w")
                    output.write("time,a,e,inc,Omega,omega,M,RESO\n")
                else:
                    output = open(histRESO_folder +
                                  "/pa_{0:d}.csv".format(i), "a")
                    output.write("{time:d},{a:.6f},{e:.6f},{I:.6f},{O:.6f},{o:.6f},{M:.6f},{RESO:d}\n"
                                 .format(time=int(output_t), a=pa['a'], e=pa['e'], I=pa['inc'], O=pa['Omega'], o=pa['omega'], M=pa['M'], RESO=self.RESOlist[i]))
                output.close()


main = Filter(5, 4, 100000, 10000000, 365000000 *
              10, 365000000*2, 182, "temp", False)
main.readFirstSnapshots()
main.loopToWindowA()
main.loopToEnd()
