%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 1.27700 7.97600 26.81300
C  1.70500 10.25900 29.27500
C  1.82700 5.41400 29.09300
C  0.85500 5.69300 24.34900
C  1.21800 10.51800 24.46400
N  1.59300 7.92000 28.89400
C  1.77500 8.96200 29.74400
C  1.81600 8.52800 31.18600
C  2.13100 6.98700 31.08800
C  1.79700 6.74800 29.57900
C  3.53500 6.42400 31.56400
C  0.41400 8.75400 31.91200
C  0.55800 8.89600 33.45600
H  1.71800 8.24400 34.22700
N  1.40400 5.81900 26.67300
C  1.63700 4.96400 27.78200
C  1.59900 3.60000 27.34800
C  1.32800 3.64600 25.93400
C  1.17100 5.11200 25.62400
C  1.72200 2.42000 28.24700
C  1.18400 2.49900 24.94000
O  1.05900 2.62000 23.68600
C  1.20700 1.06700 25.44700
N  1.13600 8.07200 24.71300
C  0.77100 6.99400 23.90500
C  0.28500 7.38700 22.48000
C  0.38900 9.00800 22.61900
C  0.90800 9.24400 23.99000
C  1.23400 6.77600 21.44100
C  -0.88200 9.73900 22.30100
C  -0.72400 10.77600 21.22600
N  1.52400 9.91600 26.73700
C  1.45400 10.87900 25.73800
C  1.65900 12.16800 26.28700
C  1.81400 12.01500 27.69100
C  1.62000 10.64100 27.94600
C  1.67300 13.48200 25.60500
C  2.15100 12.69400 28.92300
O  2.45500 13.83200 29.22300
C  2.00700 11.52600 30.01000
C  3.06400 11.59000 31.05300
O  4.23800 11.79300 30.87600
O  2.59100 11.43200 32.29900
C  3.63300 11.47700 33.36400
H  2.18200 4.70500 29.84400
H  0.85400 4.95200 23.54700
H  1.09900 11.38200 23.80700
H  2.67100 8.94100 31.72100
H  1.35600 6.44600 31.63000
H  4.07900 7.25400 32.01500
H  4.16900 5.94900 30.81700
H  3.32200 5.63000 32.28000
H  -0.18400 7.89400 31.61000
H  -0.01700 9.66700 31.50100
H  -0.39200 8.49700 33.81200
H  0.49500 9.96300 33.66800
H  0.81200 2.12400 28.76900
H  2.52200 2.67200 28.94200
H  2.18600 1.56000 27.76500
H  2.26000 0.80200 25.54100
H  0.75700 0.48100 24.64600
H  0.71900 0.92900 26.41200
H  -0.69500 6.92500 22.36000
H  1.13100 9.36700 21.90500
H  2.13800 6.31200 21.83600
H  1.41900 7.51100 20.65800
H  0.69100 5.93800 21.00300
H  -1.27500 10.30600 23.14500
H  -1.71400 9.09900 22.00700
H  -1.48700 10.70900 20.45100
H  0.24000 10.68600 20.72600
H  -0.68400 11.71900 21.77200
H  2.02800 13.48300 24.57400
H  2.51200 13.99500 26.07600
H  0.69500 13.95600 25.68900
H  1.08300 11.88900 30.46100
H  3.95200 12.49700 33.58300
H  4.44500 10.77600 33.17000
H  3.18800 11.10700 34.28800
Mg 40.61000 41.38800 26.95600
C  39.86800 43.56800 29.58900
C  41.56600 39.09100 29.09600
C  42.01100 39.61800 24.36600
C  40.07100 44.02500 24.67100
N  40.59100 41.34900 29.06300
C  40.20000 42.29100 29.98900
C  40.16900 41.71800 31.40000
C  40.93300 40.36900 31.20600
C  40.94200 40.20500 29.72600
C  42.27000 40.16100 31.88400
C  38.74600 41.61400 32.03600
C  38.70900 41.06400 33.46600
H  37.95300 41.90400 34.48500
N  41.75400 39.53700 26.76500
C  41.90400 38.67000 27.78800
C  42.49400 37.45600 27.29600
C  42.73100 37.64100 25.92300
C  42.16600 38.98800 25.63800
C  42.83400 36.26800 28.17100
C  43.30000 36.69700 24.91500
O  43.61700 36.97700 23.72800
C  43.50600 35.23000 25.27200
N  41.00200 41.78000 24.87100
C  41.64700 40.91600 24.05300
C  41.55500 41.44600 22.52200
C  40.51800 42.61700 22.64600
C  40.44500 42.80800 24.12700
C  42.99600 41.87100 22.03700
C  39.09100 42.15200 22.17700
C  38.49200 40.88200 22.76100
N  40.00500 43.44200 27.05400
C  39.84200 44.33300 26.02900
C  39.43100 45.53300 26.59500
C  39.42300 45.36800 28.01000
C  39.78400 44.06700 28.24300
C  39.15600 46.83600 25.93100
C  39.20700 45.98800 29.33700
O  38.93200 47.12500 29.70300
C  39.39600 44.78700 30.41700
C  40.24300 45.23600 31.55400
O  41.22600 46.02900 31.45500
O  39.66700 44.79200 32.66300
C  40.07000 45.50000 33.86300
H  41.79000 38.25100 29.75700
H  42.42000 39.05300 23.52500
H  39.79000 44.77400 23.92900
H  40.67200 42.47500 32.00200
H  40.36400 39.50500 31.55100
H  42.49900 40.80900 32.73000
H  43.06200 40.23200 31.13800
H  42.35700 39.09400 32.09000
H  38.09900 41.02000 31.39000
H  38.28900 42.60300 32.04400
H  39.66400 41.20800 33.97100
H  38.39800 40.01900 33.48100
H  42.71500 36.43800 29.24100
H  43.88400 35.97700 28.18900
H  42.26200 35.37200 27.92800
H  44.31400 35.17500 26.00100
H  43.97700 34.67500 24.46100
H  42.58300 34.72100 25.55200
H  41.23200 40.63200 21.87300
H  40.99300 43.46400 22.15100
H  43.79600 41.51800 22.68800
H  43.14700 42.94900 21.97800
H  43.27300 41.44900 21.07200
H  39.29200 42.04500 21.11100
H  38.30900 42.90300 22.28500
H  39.25300 40.32600 23.30900
H  37.94200 40.24300 22.07100
H  37.70100 41.27300 23.40100
H  38.42700 46.78300 25.12200
H  40.10600 47.12900 25.48500
H  38.86200 47.61200 26.63700
H  38.38700 44.62200 30.79400
H  39.19900 45.70100 34.48800
H  40.53400 46.47500 33.71600
H  40.74300 44.82600 34.39400


