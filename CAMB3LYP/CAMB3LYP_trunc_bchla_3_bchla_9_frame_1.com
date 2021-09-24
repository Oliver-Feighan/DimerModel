%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 1.39700 7.67200 26.33400
C  1.72400 9.74100 29.05100
C  2.01300 4.99100 28.38100
C  1.18200 5.57500 23.58000
C  1.22200 10.43500 24.24700
N  1.88100 7.36900 28.50400
C  1.86000 8.39700 29.40300
C  2.08800 7.91400 30.75300
C  2.40000 6.43300 30.60700
C  2.12200 6.20500 29.05700
C  3.81300 5.91300 31.04400
C  0.86900 8.21500 31.69000
C  1.18400 8.50600 33.19800
H  2.44800 7.85400 33.79200
N  1.46600 5.54900 26.07500
C  1.71900 4.62900 27.07800
C  1.70700 3.27200 26.49500
C  1.45400 3.42700 25.06400
C  1.33100 4.87900 24.84700
C  1.96100 1.99300 27.27600
C  1.22600 2.33500 24.04500
O  0.94200 2.65000 22.87600
C  1.34900 0.88300 24.32000
N  1.30400 7.94800 24.20100
C  1.17000 6.94200 23.33200
C  0.89900 7.47700 21.91600
C  0.64800 9.06900 22.12400
C  1.00600 9.14600 23.65000
C  2.06300 7.12000 20.94100
C  -0.80600 9.55500 21.87400
C  -1.00500 10.45600 20.64600
N  1.52400 9.70800 26.54100
C  1.45100 10.73400 25.59600
C  1.56500 12.05000 26.26300
C  1.69000 11.65900 27.61700
C  1.66900 10.24300 27.72900
C  1.49300 13.37100 25.67200
C  1.73800 12.19000 28.99500
O  1.68200 13.31300 29.42100
C  1.80500 10.95700 29.96500
C  2.97400 11.12000 30.84400
O  4.12200 10.98100 30.39200
O  2.68100 11.28600 32.16300
C  3.77100 11.25200 33.13300
H  2.21100 4.15900 29.06000
H  1.15400 4.90100 22.72100
H  1.03400 11.27000 23.56900
H  2.96600 8.42700 31.14700
H  1.62000 5.81700 31.05600
H  4.56300 6.57700 31.47500
H  4.27400 5.43500 30.18000
H  3.77100 5.07400 31.73900
H  0.12500 7.42500 31.58700
H  0.27300 9.01100 31.24400
H  0.33700 8.21300 33.81900
H  1.41300 9.55400 33.38900
H  2.02800 2.15000 28.35200
H  2.91300 1.54300 26.99300
H  1.08700 1.34400 27.22400
H  1.03600 0.37200 23.41000
H  0.74800 0.54700 25.16600
H  2.36800 0.65900 24.63600
H  -0.06900 7.08000 21.61000
H  1.35900 9.68000 21.56800
H  1.72500 6.67300 20.00600
H  2.74700 6.40500 21.39900
H  2.59200 8.06400 20.80800
H  -1.22600 10.10700 22.71400
H  -1.49800 8.72100 21.75500
H  -0.24900 10.22000 19.89700
H  -0.86300 11.48900 20.96200
H  -2.04500 10.34200 20.34300
H  1.26400 14.14800 26.40100
H  0.76200 13.44400 24.86700
H  2.50500 13.65900 25.38700
H  0.85200 11.05500 30.48600
H  4.47400 12.08100 33.04700
H  4.27000 10.28300 33.12600
H  3.36500 11.36600 34.13800
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


