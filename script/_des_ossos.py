import orbitplotkit as opk
from astropy.table import Table
import pandas as pd

aN = 30.0687

data_OSSOS = Table.read("Ensemble.CDS", format="ascii")
df_OSSOS = data_OSSOS.to_pandas()
df_OSSOS['q'] = df_OSSOS['a']*(1-df_OSSOS['e'])

df_OSSOS_detached = df_OSSOS[(df_OSSOS['a'] > 48) & ((df_OSSOS['cl'] == 'det') | (df_OSSOS['cl'] == 'cla'))]
# print(df_OSSOS_detached)


namelist= []
for i, MPC in enumerate(df_OSSOS_detached['MPC']):
    name = opk.Compact2MPC(MPC)
    if name != '-':
        namelist.append(name)
    else:
        namelist.append(df_OSSOS_detached.iloc[i]['object'])

df_OSSOS_detached['Name'] = namelist
df1 = df_OSSOS_detached[['Name', 'a', 'q', 'i']]
df1['Source'] = 'OSSOS'
print(len(df1.index))

# opk.fits2csv("y6_paper_table.fits", "DES_y6_paper.csv")
df_DES = pd.read_csv("DES_y6_paper.csv")
df_DES_detached = df_DES[(df_DES['a'] > 48) & ((df_DES['Class'] == 'Detached') | (df_DES['Class'] == 'Classical'))]
df_DES_detached['Name'] = df_DES_detached['MPC']
df_DES_detached['Source'] = 'DES'
df2 = df_DES_detached[['Name', 'a', 'q', 'i']]
df2['Source'] = 'DES'
print(len(df2.index))

df_det = pd.concat([df1,df2]).drop_duplicates('Name').sort_values(by = ['Name']).reset_index(drop=True)
df_det.to_csv("_DETACHED_ALL.csv")
print(df_det)