%nproc=24
%mem=175gb
#p wB97XD/Def2SVP td=(nstates=5)

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
Mg 6.73600 57.36100 41.87100
C  5.57600 54.03200 41.84400
C  9.54700 56.26000 40.29200
C  7.50300 60.62700 41.21800
C  3.48900 58.25400 43.03700
N  7.42000 55.35200 41.01200
C  6.89900 54.10900 41.23100
C  7.68800 53.00900 40.50000
C  9.08300 53.80000 40.39600
C  8.66300 55.23700 40.58100
C  10.16800 53.39100 41.42800
C  7.12500 52.63700 39.05500
C  7.43400 51.26500 38.47600
H  8.38100 51.07500 37.27200
N  8.34400 58.41300 40.90600
C  9.44600 57.74300 40.46000
C  10.47100 58.62900 40.12400
C  9.93400 59.94200 40.40300
C  8.52300 59.74000 40.83800
C  11.84700 58.22300 39.61000
C  10.74100 61.20700 40.21700
O  11.93000 61.13100 39.81300
C  10.22200 62.63600 40.60600
N  5.63100 59.10700 42.11600
C  6.19700 60.37100 41.79700
C  5.36400 61.47100 42.39400
C  3.99000 60.81500 42.65100
C  4.37000 59.30600 42.68800
C  5.80500 62.13900 43.71600
C  2.94100 61.09300 41.45000
C  2.15100 62.39300 41.46200
N  4.97100 56.31000 42.51500
C  3.75500 56.80200 42.96500
C  2.84300 55.78700 43.18500
C  3.53400 54.62400 42.80400
C  4.81800 55.00700 42.38100
C  1.45400 55.97300 43.70900
C  3.46600 53.25900 42.78000
O  2.61500 52.50900 43.21600
C  4.68000 52.84700 41.92100
C  4.09400 52.40300 40.61400
O  3.41000 53.09000 39.92100
O  4.40200 51.06200 40.37800
C  3.73000 50.54700 39.23600
H  10.54100 55.93300 39.97900
H  7.77500 61.68500 41.22200
H  2.57200 58.59300 43.52300
H  7.89200 52.14900 41.13700
H  9.49500 53.56500 39.41500
H  10.79300 52.70800 40.85200
H  9.75000 52.84400 42.27300
H  10.76800 54.24600 41.74000
H  7.52400 53.32000 38.30500
H  6.05400 52.83900 39.08200
H  6.45300 50.86600 38.21400
H  7.81100 50.65000 39.29300
H  12.72600 58.79500 39.90900
H  11.62200 58.24700 38.54400
H  12.11100 57.23800 39.99500
H  10.03600 62.67400 41.67900
H  9.32800 62.73100 39.99000
H  10.83200 63.47600 40.27400
H  5.23800 62.21700 41.61000
H  3.54200 61.14600 43.58800
H  5.15700 61.94700 44.57100
H  6.23700 63.12400 43.53800
H  6.70300 61.58800 43.99500
H  2.28600 60.23100 41.57500
H  3.43700 61.03400 40.48100
H  2.26900 62.89100 40.50000
H  2.53500 63.03900 42.25200
H  1.08300 62.20800 41.58000
H  0.80400 55.41700 43.03400
H  1.12000 57.01000 43.65700
H  1.32300 55.57400 44.71500
H  5.17600 52.01400 42.41800
H  2.78900 50.09000 39.54100
H  4.23500 49.69300 38.78300
H  3.52200 51.31600 38.49400


