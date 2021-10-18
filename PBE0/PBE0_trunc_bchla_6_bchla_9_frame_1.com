%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 17.05400 -2.28800 27.75000
C  16.54800 -0.35400 30.72600
C  18.85200 -4.52900 29.65700
C  18.01300 -3.67000 24.94900
C  15.39200 0.22200 26.03300
N  17.88600 -2.27700 29.96200
C  17.41200 -1.43700 30.93300
C  17.78200 -2.00100 32.33100
C  18.58500 -3.32700 31.99600
C  18.48400 -3.42400 30.43000
C  20.10200 -3.41100 32.45300
C  16.57800 -2.28600 33.26800
C  16.68400 -1.85000 34.68600
H  18.05300 -1.52300 35.39600
N  18.30500 -3.91800 27.31000
C  18.85900 -4.81900 28.21900
C  19.64100 -5.82500 27.47500
C  19.29300 -5.63200 26.08300
C  18.51600 -4.37100 26.05800
C  20.53100 -6.87100 28.07500
C  19.63900 -6.46300 24.82500
O  19.04500 -6.38400 23.72900
C  20.63000 -7.54600 24.97700
N  16.79100 -1.75400 25.84400
C  17.37300 -2.49900 24.84200
C  16.84500 -2.07000 23.44500
C  15.70300 -1.04100 23.86200
C  15.95300 -0.84300 25.34700
C  17.88000 -1.47900 22.40300
C  14.26500 -1.46100 23.60200
C  13.50200 -0.61900 22.55800
N  15.98100 -0.53800 28.25500
C  15.38600 0.39400 27.43800
C  14.87800 1.51400 28.24100
C  15.27300 1.20700 29.53100
C  15.93800 -0.03900 29.49700
C  14.17200 2.75100 27.79300
C  15.18900 1.69900 30.85900
O  14.57400 2.68200 31.36200
C  15.98100 0.75800 31.66400
C  17.01200 1.50300 32.46800
O  17.97500 2.03100 32.02700
O  16.72000 1.42300 33.76600
C  17.36800 2.29700 34.68400
H  19.44500 -5.28100 30.18200
H  18.21500 -4.04600 23.94300
H  14.94000 1.04300 25.47400
H  18.42200 -1.24900 32.79400
H  18.07000 -4.18900 32.41900
H  20.43900 -2.60400 33.10300
H  20.70600 -3.31100 31.55000
H  20.31100 -4.38200 32.90300
H  16.34900 -3.34900 33.35200
H  15.66000 -1.80400 32.93300
H  16.30900 -2.63700 35.34100
H  16.04900 -0.97200 34.80100
H  21.52400 -6.79800 27.63200
H  20.09600 -7.83300 27.80400
H  20.47400 -6.87000 29.16400
H  21.57900 -7.08500 25.25100
H  20.75500 -7.96300 23.97800
H  20.26200 -8.35300 25.61000
H  16.37000 -2.96000 23.03100
H  15.88800 -0.12500 23.30100
H  18.88000 -1.46500 22.83500
H  17.57100 -0.51800 21.99000
H  17.86400 -2.13300 21.53100
H  13.63500 -1.61600 24.47800
H  14.30700 -2.45900 23.16700
H  12.47900 -0.33500 22.80400
H  13.29800 -1.28100 21.71600
H  13.99200 0.30400 22.25100
H  13.75500 2.65600 26.79000
H  14.88900 3.56000 27.65400
H  13.49700 3.03000 28.60200
H  15.27200 0.27100 32.33300
H  16.91900 3.28000 34.53500
H  18.41700 2.43200 34.42000
H  17.17400 2.00000 35.71500
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


