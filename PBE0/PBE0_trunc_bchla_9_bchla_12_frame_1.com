%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 35.71600 1.25500 29.88400
C  33.26000 2.26900 32.06200
C  37.97500 1.18400 32.63200
C  38.34600 0.78200 27.81300
C  33.67000 1.64600 27.25500
N  35.61600 1.68700 32.10200
C  34.47500 2.10300 32.73400
C  34.68900 2.20300 34.24600
C  36.27600 2.03500 34.33300
C  36.65900 1.61200 32.94500
C  37.02900 3.30200 34.88200
C  33.91200 1.08400 35.07900
C  33.16700 1.39000 36.35600
H  33.84300 2.22100 37.38800
N  37.88700 0.95800 30.15700
C  38.55000 0.78800 31.32200
C  39.89100 0.36100 31.10100
C  40.02700 0.20000 29.68600
C  38.73000 0.72700 29.16600
C  40.80200 0.22400 32.27800
C  41.26400 -0.16900 28.88500
O  41.33400 -0.00500 27.66300
C  42.48400 -0.74400 29.58300
N  36.02000 1.36700 27.77900
C  37.17800 1.00100 27.16900
C  36.98100 0.89100 25.65800
C  35.41500 1.00100 25.51900
C  34.99600 1.35200 26.93300
C  37.78500 1.78800 24.67100
C  34.66500 -0.24600 24.88900
C  34.10400 -1.33400 25.78000
N  33.89900 1.92200 29.64700
C  33.07800 1.91600 28.54800
C  31.77800 2.41500 28.81300
C  31.77700 2.52500 30.20000
C  33.07800 2.19300 30.65800
C  30.56200 2.67700 27.93400
C  31.02900 2.87500 31.32700
O  29.88700 3.29400 31.37200
C  31.95500 2.75800 32.66400
C  32.02000 4.07800 33.27900
O  32.54800 5.09200 32.78100
O  31.38700 4.14700 34.48000
C  31.10000 5.46600 35.05500
H  38.54900 1.11500 33.55800
H  39.09900 0.35600 27.14700
H  32.90400 1.43100 26.50800
H  34.28100 3.15500 34.58800
H  36.55500 1.25900 35.04600
H  37.56800 3.00200 35.78000
H  36.24700 4.05100 35.00800
H  37.83800 3.67300 34.25200
H  34.72200 0.38900 35.30100
H  33.30700 0.45900 34.42300
H  32.77000 0.47000 36.78500
H  32.33500 2.02700 36.05800
H  41.31400 -0.72500 32.12300
H  40.23000 0.11700 33.20000
H  41.52900 1.02400 32.42200
H  42.25000 -1.48700 30.34600
H  43.15800 0.06500 29.86800
H  42.99900 -1.30000 28.79900
H  37.38700 -0.09500 25.43100
H  35.10800 1.75000 24.79000
H  37.12800 2.52600 24.21100
H  38.16400 1.12200 23.89600
H  38.57800 2.24100 25.26500
H  35.37100 -0.71500 24.20500
H  33.81300 0.09600 24.30100
H  34.44300 -2.32400 25.47400
H  33.01400 -1.34700 25.82000
H  34.48200 -1.09900 26.77500
H  30.79100 2.37600 26.91100
H  30.30400 3.73500 27.98700
H  29.72700 2.09100 28.31700
H  31.32600 2.10900 33.27400
H  30.27200 5.47100 35.76400
H  31.01600 6.28300 34.33900
H  31.93600 5.66400 35.72600
Mg 47.67100 15.89000 28.45400
C  45.34400 15.65500 31.04800
C  49.48000 18.07500 30.51500
C  49.68700 16.65900 25.89600
C  45.30600 14.62100 26.30000
N  47.50700 16.71700 30.64000
C  46.40900 16.47700 31.40500
C  46.52000 17.23500 32.78700
C  47.75800 18.22000 32.48000
C  48.27100 17.69900 31.12700
C  47.16600 19.65600 32.25100
C  46.90400 16.16800 33.87700
C  47.20400 16.68600 35.31000
H  48.31600 16.29900 36.23900
N  49.43000 17.01400 28.32000
C  50.06900 17.82700 29.27100
C  51.20000 18.39100 28.66300
C  51.38300 17.85100 27.39700
C  50.14500 17.12800 27.15700
C  52.11500 19.27100 29.45000
C  52.51500 18.07600 26.34100
O  52.57400 17.52300 25.24800
C  53.78800 18.86100 26.66900
N  47.44200 15.76700 26.39600
C  48.48400 16.01700 25.57200
C  48.26000 15.53000 24.12700
C  46.80100 15.00400 24.22400
C  46.51000 15.00800 25.71400
C  48.61100 16.58300 23.02400
C  46.48600 13.65000 23.53100
C  46.96600 12.43800 24.21000
N  45.75600 15.16700 28.59900
C  44.93800 14.69700 27.68600
C  43.69800 14.18000 28.27400
C  43.86300 14.56100 29.64400
C  45.08900 15.18100 29.78600
C  42.77800 13.18000 27.74800
C  43.21100 14.49100 30.94000
O  42.19200 13.93500 31.28200
C  44.14500 15.27600 31.86900
C  43.41900 16.43800 32.55500
O  43.17000 17.49000 32.02600
O  43.04200 16.05500 33.82000
C  42.33900 17.13300 34.55300
H  50.09800 18.64700 31.20900
H  50.41400 16.70200 25.08200
H  44.81200 14.14000 25.45400
H  45.63400 17.83000 33.01000
H  48.52200 18.24600 33.25800
H  47.13800 20.01300 31.22200
H  47.84000 20.42900 32.62100
H  46.16600 19.79300 32.66200
H  47.75800 15.61000 33.49400
H  46.09800 15.43400 33.89800
H  46.24800 16.58700 35.82500
H  47.42100 17.73200 35.09700
H  51.80300 19.48300 30.47300
H  52.34700 20.14300 28.83700
H  53.02500 18.67700 29.53900
H  53.61800 19.92500 26.83400
H  54.61200 18.70300 25.97300
H  54.19200 18.41700 27.57900
H  49.00300 14.73300 24.13500
H  46.20400 15.80500 23.78800
H  49.40000 17.25400 23.36200
H  47.72600 17.18100 22.80500
H  48.98600 16.12600 22.10800
H  46.99700 13.77600 22.57700
H  45.42300 13.59700 23.29500
H  46.25100 12.14600 24.97900
H  47.97200 12.58100 24.60400
H  46.98800 11.57800 23.54000
H  43.37600 12.59600 27.04800
H  41.93300 13.49000 27.13300
H  42.41600 12.61300 28.60500
H  44.52600 14.58100 32.61800
H  41.57800 16.67700 35.18600
H  41.70900 17.79000 33.95300
H  43.02300 17.66800 35.21100


