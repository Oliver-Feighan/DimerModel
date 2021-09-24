%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 2.79300 -0.99900 44.50000
C  6.05500 -0.17400 43.38400
C  1.60600 1.42400 42.43900
C  -0.30000 -2.33900 44.94400
C  4.15400 -3.74300 46.12200
N  3.83100 0.31200 42.98000
C  5.15200 0.58300 42.74000
C  5.32600 1.71500 41.79700
C  3.88200 2.40600 41.78900
C  3.05100 1.33800 42.41400
C  3.87900 3.72200 42.55000
C  5.83700 1.13000 40.47100
C  7.17200 1.73800 39.86600
H  7.22000 2.30200 38.41400
N  0.88600 -0.54700 43.86000
C  0.60800 0.57500 43.06700
C  -0.81500 0.72700 43.01300
C  -1.42800 -0.37100 43.70600
C  -0.26400 -1.10400 44.27200
C  -1.44600 1.92200 42.31000
C  -2.89700 -0.69100 43.87000
O  -3.66200 0.20700 43.41300
C  -3.55500 -1.99100 44.41300
N  1.99700 -2.80800 45.42300
C  0.68200 -3.10900 45.50100
C  0.40600 -4.46500 46.20300
C  1.84700 -4.85700 46.68000
C  2.74500 -3.73600 45.97300
C  -0.59100 -4.26300 47.33000
C  2.30300 -6.22900 46.28700
C  2.45600 -7.17500 47.44700
N  4.64000 -1.73900 44.89100
C  5.11000 -2.84400 45.61200
C  6.47200 -2.80600 45.78400
C  6.89300 -1.74700 44.95400
C  5.75400 -1.18800 44.37900
C  7.21700 -3.66700 46.78200
C  8.08600 -0.96000 44.54300
O  9.22400 -1.08600 44.87400
C  7.50800 0.01200 43.44600
C  8.17400 1.35500 43.65200
O  9.07700 1.80400 42.94600
O  7.57500 2.04000 44.73700
C  8.09900 3.39600 44.89800
H  1.22700 2.33000 41.96100
H  -1.28000 -2.75100 45.19600
H  4.63400 -4.67900 46.41300
H  6.07700 2.40100 42.18900
H  3.61200 2.61400 40.75400
H  4.50800 3.70600 43.43900
H  2.87200 3.89300 42.93200
H  4.11600 4.58800 41.93200
H  5.03800 1.25300 39.74000
H  5.94100 0.04500 40.47900
H  7.81400 0.85800 39.91300
H  7.65200 2.40300 40.58300
H  -1.04200 2.04400 41.30500
H  -1.16400 2.79800 42.89400
H  -2.51800 1.88200 42.11200
H  -3.31300 -1.96500 45.47500
H  -3.08800 -2.84600 43.92400
H  -4.62700 -2.03100 44.21800
H  -0.01900 -5.22600 45.54900
H  1.95500 -4.64800 47.74400
H  -1.44700 -3.72600 46.92100
H  -0.11100 -3.66200 48.10300
H  -0.92100 -5.12200 47.91300
H  3.28200 -6.24400 45.80800
H  1.69000 -6.76000 45.55900
H  1.86400 -8.06400 47.23300
H  2.03500 -6.78100 48.37200
H  3.49000 -7.42400 47.68800
H  8.18600 -3.18800 46.92600
H  7.36400 -4.68500 46.42200
H  6.59100 -3.65900 47.67400
H  7.98500 -0.36600 42.54100
H  7.35600 4.03900 45.37000
H  8.23900 3.88500 43.93500
H  9.04800 3.42700 45.43400
Mg 40.80400 41.48400 27.15800
C  39.93800 43.96600 29.51800
C  41.67100 39.43500 29.59200
C  42.41200 39.59100 24.81200
C  40.80000 44.20800 24.65700
N  40.84000 41.74700 29.31900
C  40.23400 42.70800 30.04600
C  40.11400 42.23500 31.50800
C  40.62500 40.83000 31.52200
C  41.12800 40.63200 30.06000
C  41.62900 40.55000 32.71400
C  38.74000 42.40700 32.23800
C  38.20000 41.35400 33.29800
H  37.58800 41.87900 34.60500
N  41.67300 39.62800 27.19200
C  41.94100 38.92200 28.34300
C  42.50600 37.66000 27.95800
C  42.79100 37.65600 26.53700
C  42.34000 39.00600 26.13300
C  42.84100 36.69600 29.02300
C  43.34600 36.53900 25.61900
O  43.54900 36.75600 24.47000
C  43.39500 35.12400 26.09600
N  41.54800 41.90500 25.04000
C  42.04300 40.86300 24.33500
C  42.34300 41.24800 22.85100
C  41.79900 42.70700 22.76800
C  41.36100 43.00200 24.22400
C  43.78600 41.15300 22.43200
C  40.63300 43.00600 21.76200
C  39.24100 42.61000 22.12200
N  40.45200 43.63600 26.97900
C  40.44000 44.51600 25.93500
C  39.94600 45.79200 26.39200
C  39.72700 45.66500 27.78900
C  40.09700 44.34100 28.10100
C  39.86300 47.06500 25.62600
C  39.33400 46.33800 28.99400
O  38.92600 47.46900 29.23200
C  39.40400 45.25700 30.16400
C  40.16300 45.75100 31.38400
O  41.33900 46.00900 31.38400
O  39.30000 45.80400 32.46700
C  39.93400 46.30500 33.68700
H  41.79100 38.62500 30.31500
H  42.89000 38.90600 24.10900
H  40.69600 44.98800 23.90000
H  40.78700 42.96300 31.96100
H  39.78000 40.17800 31.74500
H  41.64600 41.41100 33.38200
H  42.57700 40.23500 32.27800
H  41.25100 39.76900 33.37400
H  37.98000 42.45900 31.45800
H  38.66600 43.37100 32.74200
H  39.06500 40.72800 33.51800
H  37.46000 40.71500 32.81600
H  43.75800 36.20700 28.69500
H  42.05100 35.94600 29.03300
H  42.98100 37.00500 30.05900
H  44.39900 35.05700 26.51300
H  43.13900 34.46100 25.26900
H  42.74200 34.89700 26.93900
H  41.71500 40.69400 22.15300
H  42.49200 43.53200 22.60600
H  44.38600 40.63600 23.18200
H  44.24800 42.13400 22.31800
H  43.78100 40.59500 21.49600
H  40.84900 42.53700 20.80200
H  40.70100 44.05800 21.48400
H  39.26600 41.80800 22.86000
H  38.56300 42.19300 21.37700
H  38.76600 43.53300 22.45600
H  40.81300 47.42600 25.23300
H  39.70300 47.84000 26.37600
H  39.09100 47.02600 24.85700
H  38.36800 45.04000 30.42500
H  39.84500 45.31900 34.14200
H  39.21300 47.03100 34.06200
H  40.94500 46.67900 33.52400


