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
Mg 29.51000 58.90800 41.21400
C  26.38900 57.48600 40.17300
C  31.07800 56.05600 39.82600
C  32.44700 60.51200 41.46200
C  27.84100 61.85400 41.88400
N  28.82700 56.95900 40.12100
C  27.51800 56.70000 39.82500
C  27.43600 55.53200 38.85500
C  28.91100 54.97600 39.03000
C  29.71500 56.02800 39.69400
C  29.08100 53.62400 39.80100
C  26.96300 56.03100 37.45200
C  25.85600 55.22200 36.65900
H  26.19700 54.60200 35.29200
N  31.43900 58.34900 40.91600
C  31.91100 57.13400 40.35300
C  33.36200 57.11100 40.33800
C  33.74700 58.43000 40.81400
C  32.51100 59.13800 41.09000
C  34.11400 55.98100 39.76900
C  35.21500 58.86100 40.97700
O  36.08800 58.03800 40.57400
C  35.69500 60.19600 41.58500
N  30.14600 61.00100 41.36300
C  31.38600 61.37000 41.62100
C  31.40400 62.84700 41.84900
C  29.90900 63.22000 42.03600
C  29.21400 61.95000 41.65200
C  32.31600 63.46200 42.99500
C  29.40100 64.48800 41.17400
C  28.93400 65.75600 41.93000
N  27.56300 59.51200 41.26600
C  27.03200 60.71300 41.55100
C  25.61900 60.69600 41.44600
C  25.33000 59.42300 40.94900
C  26.56000 58.73200 40.80800
C  24.61100 61.72700 41.76100
C  24.21900 58.63400 40.51400
O  23.05400 58.95300 40.43300
C  24.92200 57.26000 39.91000
C  24.42400 56.08500 40.63500
O  24.27500 56.00600 41.86100
O  24.19400 55.10700 39.69100
C  23.68300 53.83100 40.18900
H  31.62000 55.20100 39.41900
H  33.35900 61.11200 41.46000
H  27.25300 62.67600 42.29800
H  26.78700 54.79200 39.32400
H  29.19500 54.98000 37.97700
H  29.75800 53.81000 40.63600
H  29.53300 52.82700 39.21100
H  28.17800 53.23200 40.27000
H  27.91200 55.97700 36.91900
H  26.65200 57.07400 37.39700
H  25.02500 55.92000 36.55900
H  25.42300 54.48900 37.33900
H  33.71100 55.81900 38.76900
H  33.89700 55.07100 40.32800
H  35.20300 56.02600 39.79800
H  35.47200 60.98900 40.87100
H  36.68900 59.96900 41.96900
H  35.05300 60.30300 42.46000
H  31.80800 63.29500 40.94200
H  29.82800 63.26800 43.12200
H  32.55800 64.46300 42.63700
H  33.18900 62.82200 43.12200
H  31.69200 63.52000 43.88700
H  28.59500 64.07600 40.56700
H  30.17400 64.71200 40.43900
H  27.96900 66.12300 41.58300
H  29.69100 66.53100 41.80700
H  28.89100 65.63000 43.01200
H  23.61700 61.29700 41.88000
H  24.71500 62.31200 40.84700
H  25.04500 62.40800 42.49300
H  24.83500 57.21500 38.82500
H  22.80800 53.96700 40.82500
H  24.47100 53.30400 40.72900
H  23.37200 53.22500 39.33800


