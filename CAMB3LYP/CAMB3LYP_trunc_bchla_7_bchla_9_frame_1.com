%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 26.07000 0.32300 29.17500
C  27.93000 0.00300 32.17200
C  23.23700 0.37300 31.22200
C  24.20000 0.73500 26.50800
C  28.89900 -0.55900 27.35600
N  25.68800 0.23300 31.42600
C  26.58700 0.25300 32.44900
C  25.99300 0.25100 33.86800
C  24.43500 0.24800 33.45400
C  24.48700 0.27000 31.95500
C  23.63400 -1.02300 33.94300
C  26.39900 1.50900 34.70300
C  27.30600 1.47800 35.92700
H  27.07500 2.41400 37.10600
N  24.00400 0.60600 28.87400
C  22.98600 0.52400 29.83500
C  21.72400 0.61300 29.13500
C  21.97300 0.88000 27.79100
C  23.44300 0.80400 27.68300
C  20.46800 0.64800 29.93700
C  21.02800 1.08200 26.60700
O  21.39200 1.28100 25.44400
C  19.58500 1.18400 26.91000
N  26.42900 0.02500 27.20000
C  25.57000 0.46200 26.23000
C  26.30100 0.60400 24.85200
C  27.70600 0.10100 25.15700
C  27.68000 -0.24400 26.67200
C  25.55900 -0.10800 23.67400
C  28.73300 1.18300 24.79000
C  29.96300 0.77300 24.06500
N  28.00500 -0.28200 29.61000
C  29.05600 -0.58000 28.78000
C  30.24500 -0.86800 29.52500
C  29.92400 -0.56000 30.93700
C  28.52000 -0.19800 30.89800
C  31.58600 -1.33700 28.93500
C  30.40700 -0.48700 32.33900
O  31.55800 -0.58600 32.79700
C  29.11500 -0.11100 33.17500
C  29.05800 -1.28400 34.10300
O  28.53100 -2.37500 33.92200
O  29.56000 -0.97800 35.26800
C  29.75400 -2.03800 36.28700
H  22.46700 0.52900 31.98000
H  23.69600 0.81600 25.54300
H  29.79300 -0.76400 26.76400
H  26.26300 -0.69600 34.33700
H  23.82700 1.11200 33.72400
H  24.21800 -1.88400 34.26700
H  22.99700 -1.39000 33.13900
H  22.99100 -0.72400 34.77100
H  25.46900 1.94600 35.06700
H  26.90100 2.23700 34.06700
H  28.23400 1.84700 35.49000
H  27.27900 0.42100 36.19200
H  19.92700 1.58400 29.80500
H  20.67800 0.53600 31.00100
H  19.76300 -0.13800 29.66800
H  19.43300 1.90200 27.71600
H  19.18900 0.27500 27.36300
H  18.97400 1.53500 26.07900
H  26.25300 1.68100 24.69200
H  27.96600 -0.74800 24.52500
H  24.56100 -0.46700 23.92500
H  26.08600 -0.97800 23.28400
H  25.56500 0.53500 22.79400
H  29.00400 1.58600 25.76500
H  28.34900 1.97700 24.15000
H  30.78000 1.14300 24.68400
H  29.93300 1.13400 23.03700
H  30.14200 -0.30200 24.03500
H  31.60000 -2.42600 28.91000
H  32.38900 -0.93100 29.55000
H  31.71700 -1.01400 27.90200
H  29.35100 0.80400 33.71800
H  28.83300 -2.53300 36.59300
H  30.27100 -1.70300 37.18600
H  30.41000 -2.86500 36.01700
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


