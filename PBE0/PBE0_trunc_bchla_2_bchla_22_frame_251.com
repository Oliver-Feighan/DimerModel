%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 3.63500 0.41800 44.38500
C  6.38200 2.41600 44.15700
C  1.85600 2.91600 42.76500
C  1.03400 -1.70200 44.27100
C  5.65900 -2.27700 45.39800
N  4.15300 2.34100 43.32400
C  5.35700 2.93700 43.43600
C  5.44100 4.20100 42.66100
C  3.82600 4.53100 42.53700
C  3.19500 3.17000 42.89900
C  3.38200 5.72400 43.43400
C  6.08500 3.93300 41.31000
C  5.45900 2.81500 40.35400
H  5.36900 3.04800 38.85100
N  1.72200 0.57100 43.63100
C  1.11700 1.74400 43.11300
C  -0.27500 1.47100 42.97900
C  -0.55500 0.16000 43.34400
C  0.72000 -0.40500 43.81700
C  -1.17200 2.58100 42.47000
C  -1.87000 -0.57100 43.21500
O  -2.81800 0.09700 42.73200
C  -2.07300 -1.97000 43.39000
N  3.45400 -1.75300 44.57800
C  2.23500 -2.36400 44.54600
C  2.35900 -3.87900 44.83600
C  3.87500 -4.03200 45.19700
C  4.39100 -2.62000 44.97700
C  1.41200 -4.36800 45.91100
C  4.59000 -5.14900 44.26500
C  5.56000 -6.06500 45.02000
N  5.66900 0.04500 44.75000
C  6.32000 -0.99300 45.33600
C  7.63200 -0.52000 45.68500
C  7.75100 0.85800 45.29200
C  6.52700 1.13100 44.70400
C  8.77000 -1.38000 46.11700
C  8.57600 2.01500 45.15500
O  9.73800 2.10000 45.57200
C  7.71800 3.05600 44.46000
C  7.61900 4.24800 45.33700
O  8.63700 4.80700 45.69500
O  6.34700 4.66100 45.59200
C  6.23700 5.85100 46.43800
H  1.24500 3.78600 42.51600
H  0.13600 -2.29700 44.44900
H  6.31200 -3.00800 45.88000
H  6.00100 5.03200 43.09000
H  3.58600 4.93100 41.55100
H  2.93500 6.50500 42.82000
H  4.28700 6.13300 43.88400
H  2.70200 5.45700 44.24400
H  7.14700 3.76700 41.49200
H  6.12500 4.93200 40.87600
H  4.43400 2.52000 40.57900
H  6.17000 1.99900 40.48000
H  -1.10900 2.52700 41.38300
H  -0.76100 3.53900 42.78600
H  -2.21500 2.54500 42.78600
H  -2.88500 -2.34300 42.76500
H  -2.15000 -2.27300 44.43500
H  -1.20000 -2.43300 42.93200
H  2.14200 -4.33200 43.86800
H  3.97500 -4.21900 46.26600
H  0.70500 -3.61400 46.25600
H  1.77700 -4.91000 46.78300
H  0.73700 -5.01500 45.35000
H  5.15900 -4.73200 43.43400
H  3.86400 -5.76700 43.73700
H  5.09600 -7.04700 45.10800
H  5.65500 -5.69200 46.04000
H  6.52700 -6.16100 44.52600
H  8.71100 -1.83500 47.10500
H  9.65100 -0.74100 46.04600
H  8.91100 -2.23600 45.45700
H  8.25000 3.32500 43.54700
H  6.78900 6.65600 45.95200
H  6.55900 5.68400 47.46600
H  5.19100 6.09800 46.61900
Mg 8.73000 48.00700 24.72600
C  6.54700 48.35400 27.54300
C  11.38500 48.83200 26.74200
C  10.78900 48.16000 22.05300
C  5.85600 48.14800 22.68600
N  8.91100 48.51400 26.88400
C  7.92400 48.53400 27.88200
C  8.57600 48.79800 29.21500
C  9.96600 49.43600 28.79200
C  10.09000 48.94100 27.32100
C  9.83300 50.99700 28.89900
C  8.62100 47.55100 30.16700
C  8.25400 47.92900 31.51300
H  7.91700 46.69900 32.49700
N  10.84000 48.27700 24.55100
C  11.74400 48.56500 25.49500
C  13.05300 48.54500 24.91800
C  12.92300 48.37200 23.55600
C  11.42000 48.19500 23.29600
C  14.33000 48.90700 25.74100
C  14.02400 48.24600 22.48300
O  13.70400 48.13100 21.32700
C  15.48100 48.17400 22.87500
N  8.33400 48.41200 22.66100
C  9.38600 48.28200 21.76100
C  8.87500 48.15200 20.36300
C  7.28100 48.12500 20.55900
C  7.12600 48.16700 22.06100
C  9.40100 49.20400 19.44600
C  6.46300 46.89200 19.88400
C  5.44500 47.20900 18.86100
N  6.61400 48.15800 24.98900
C  5.61400 48.20500 24.04400
C  4.27700 48.23000 24.78200
C  4.59000 48.17400 26.15500
C  6.02700 48.18500 26.22900
C  2.90700 48.32700 24.09700
C  3.98700 48.26800 27.45000
O  2.78100 48.32100 27.78100
C  5.21700 48.56000 28.40400
C  5.00500 47.63900 29.64300
O  5.64000 46.57500 29.76100
O  4.15000 48.22900 30.51400
C  3.80400 47.24100 31.59300
H  12.25300 48.87100 27.40400
H  11.43600 48.18800 21.17300
H  4.90900 48.09600 22.14500
H  8.00400 49.60800 29.66700
H  10.77600 49.10700 29.44200
H  9.23900 51.44000 28.10000
H  10.83100 51.43100 28.94300
H  9.35800 51.24300 29.84900
H  9.67300 47.28200 30.26500
H  8.05700 46.70900 29.76500
H  7.29600 48.44900 31.52800
H  9.07600 48.51700 31.92200
H  14.96900 48.03700 25.59000
H  14.16900 49.02000 26.81300
H  14.74800 49.84700 25.38100
H  16.11600 48.17200 21.98900
H  15.42200 47.21200 23.38400
H  15.69700 49.01800 23.53000
H  9.23500 47.16500 20.07200
H  6.89900 49.04100 20.10900
H  10.07700 48.66800 18.78100
H  10.03800 49.93000 19.95100
H  8.64800 49.76300 18.89000
H  5.91200 46.49600 20.73700
H  7.19200 46.14700 19.56800
H  4.44900 46.91600 19.19200
H  5.77800 46.47600 18.12700
H  5.55100 48.21100 18.44300
H  2.66500 49.38900 24.07600
H  2.13000 47.78700 24.63800
H  2.97200 47.86900 23.11000
H  5.09400 49.61500 28.64500
H  3.12300 46.44700 31.28700
H  3.39400 47.68200 32.50100
H  4.75700 46.82700 31.92300


