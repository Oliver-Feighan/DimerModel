%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg -5.36100 24.76000 27.02500
C  -3.64200 26.63800 29.43600
C  -6.40000 22.65900 29.54700
C  -6.94500 22.99100 24.73300
C  -4.07500 26.91800 24.56900
N  -5.09300 24.76200 29.22400
C  -4.33400 25.60400 30.03700
C  -4.64200 25.43000 31.57600
C  -5.37000 24.10100 31.51700
C  -5.66000 23.80300 29.98800
C  -4.43700 22.92100 32.05400
C  -5.65800 26.61300 31.90200
C  -5.45100 27.38800 33.27700
H  -4.92400 26.43700 34.32100
N  -6.35800 23.02700 27.11800
C  -6.68000 22.31300 28.19300
C  -7.48400 21.16300 27.87400
C  -7.81800 21.31800 26.51400
C  -6.97100 22.47200 26.02900
C  -7.92400 20.11100 28.82500
C  -8.85200 20.60400 25.66800
O  -9.21200 21.01500 24.55900
C  -9.41800 19.25800 26.12900
N  -5.43800 24.92500 24.91100
C  -6.26300 24.05500 24.23800
C  -6.35200 24.50400 22.78100
C  -5.50900 25.87000 22.72900
C  -5.01000 25.97700 24.14600
C  -5.90500 23.33400 21.83500
C  -6.23300 27.10700 22.22100
C  -5.47900 27.85000 21.07400
N  -4.10000 26.41500 26.99800
C  -3.63700 27.19900 25.89500
C  -2.75500 28.21500 26.34700
C  -2.74300 28.01700 27.73800
C  -3.53100 26.89000 28.05100
C  -2.03600 29.28900 25.59600
C  -2.22900 28.50400 28.97800
O  -1.46700 29.40000 29.25600
C  -2.89200 27.68300 30.12000
C  -1.83900 27.24500 31.07200
O  -0.96900 26.39800 30.89000
O  -2.04600 27.91300 32.21100
C  -1.18400 27.58800 33.32800
H  -6.90800 22.00000 30.25400
H  -7.48200 22.35700 24.02300
H  -3.70300 27.64900 23.84800
H  -3.68500 25.58900 32.07200
H  -6.31300 24.16400 32.05900
H  -3.44500 23.27200 32.33900
H  -4.33700 22.17800 31.26300
H  -4.99100 22.51200 32.89900
H  -6.69800 26.29000 31.93800
H  -5.59700 27.30300 31.06000
H  -6.41300 27.65600 33.71400
H  -4.98200 28.36000 33.12200
H  -9.00400 20.22200 28.93100
H  -7.47900 20.25200 29.80900
H  -7.70000 19.10000 28.48600
H  -9.79600 18.88600 25.17700
H  -10.18200 19.34000 26.90300
H  -8.70400 18.49400 26.43600
H  -7.35300 24.74000 22.42000
H  -4.63900 25.71000 22.09200
H  -6.65200 23.29700 21.04200
H  -5.78500 22.37600 22.34100
H  -5.04100 23.62000 21.23400
H  -6.46800 27.83800 22.99400
H  -7.10900 26.70500 21.71000
H  -4.60800 27.27900 20.75200
H  -5.00600 28.76500 21.43000
H  -6.13000 28.01900 20.21600
H  -2.85300 29.97800 25.38000
H  -1.55300 28.80400 24.74800
H  -1.32200 29.82700 26.21800
H  -3.47700 28.43300 30.65400
H  -0.13900 27.68100 33.03100
H  -1.23800 26.56800 33.71000
H  -1.41700 28.31200 34.10900


