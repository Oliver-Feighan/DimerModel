%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg 35.28300 49.28100 25.05100
C  35.19500 47.66600 28.21300
C  33.70700 52.00200 26.67100
C  34.96700 50.74000 22.15000
C  36.56200 46.36300 23.65600
N  34.51200 49.78400 27.24100
C  34.71400 48.96400 28.32400
C  34.31200 49.60400 29.64200
C  33.71300 50.96400 29.08000
C  34.02200 50.95100 27.57500
C  32.25500 51.13300 29.49700
C  35.54400 49.80300 30.55200
C  35.36900 49.62200 32.10800
H  33.94300 49.72800 32.69800
N  34.52100 51.15700 24.48200
C  33.85000 52.07800 25.28700
C  33.47900 53.23000 24.45000
C  33.93900 52.90600 23.13200
C  34.54200 51.56500 23.25500
C  32.73100 54.35200 24.94000
C  33.82500 53.73400 21.85600
O  34.18700 53.27500 20.76100
C  33.44700 55.16200 21.91300
N  35.61000 48.60100 23.15100
C  35.50000 49.43800 22.11600
C  35.83700 48.71600 20.77100
C  36.52000 47.40900 21.28500
C  36.22400 47.41000 22.75200
C  34.60600 48.52300 19.81400
C  38.05900 47.61300 21.01100
C  38.99100 48.19700 22.12800
N  35.88400 47.43300 25.74900
C  36.37900 46.32800 25.03300
C  36.56500 45.20400 25.98400
C  36.06300 45.68400 27.18800
C  35.69100 47.00600 27.06500
C  36.87500 43.87200 25.49400
C  35.71300 45.23100 28.53200
O  35.83100 44.10500 29.00200
C  35.09700 46.50900 29.27300
C  33.71900 46.24000 29.75100
O  32.68000 46.50000 29.15100
O  33.78400 45.64900 31.02700
C  32.59200 45.22000 31.72300
H  33.45300 52.96300 27.12400
H  34.90900 51.17800 21.15100
H  37.00400 45.52600 23.11200
H  33.56100 49.02400 30.17900
H  34.21800 51.75900 29.62800
H  32.03300 52.05800 30.02800
H  31.83200 50.33100 30.10200
H  31.68100 51.22500 28.57400
H  35.98600 50.78700 30.39600
H  36.22100 49.06000 30.12900
H  35.89600 50.41200 32.64200
H  35.97000 48.74800 32.35900
H  33.30100 55.24000 24.66600
H  32.41200 54.18900 25.96900
H  31.80500 54.22200 24.37900
H  32.39300 55.32500 22.14000
H  33.61200 55.61800 20.93700
H  33.97000 55.60700 22.75900
H  36.50600 49.47100 20.35800
H  36.19100 46.53200 20.72700
H  34.58000 49.32300 19.07400
H  33.72600 48.44000 20.45100
H  34.76400 47.56600 19.31700
H  38.19100 48.07200 20.03100
H  38.44400 46.59300 21.03000
H  39.60900 48.91300 21.58700
H  39.65000 47.44800 22.56900
H  38.49600 48.74500 22.92900
H  36.19400 43.63700 24.67700
H  36.66800 43.18600 26.31600
H  37.94500 43.75700 25.32300
H  35.75600 46.63900 30.13200
H  32.37400 44.16000 31.59500
H  31.77700 45.93200 31.59200
H  32.87100 45.44400 32.75300


