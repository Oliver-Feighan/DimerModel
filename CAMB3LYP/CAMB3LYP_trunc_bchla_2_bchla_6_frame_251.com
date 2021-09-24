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
Mg 17.10100 -2.07700 27.77000
C  16.45500 -0.19100 30.58100
C  19.23600 -3.94700 29.68100
C  17.96100 -3.74600 24.96300
C  15.37700 0.23000 25.82200
N  17.83400 -1.91300 29.91400
C  17.31200 -1.24300 30.86800
C  17.86400 -1.55200 32.35400
C  18.99500 -2.58100 31.95700
C  18.76000 -2.80300 30.38800
C  20.50900 -2.04800 32.15900
C  16.79200 -2.13400 33.37400
C  16.96200 -1.78000 34.89400
H  18.25700 -1.14000 35.29400
N  18.31900 -3.63600 27.38700
C  19.10500 -4.38300 28.35000
C  19.74900 -5.47800 27.63200
C  19.48000 -5.36700 26.24500
C  18.51900 -4.24500 26.14900
C  20.50800 -6.46400 28.44500
C  19.91900 -6.24400 25.04900
O  19.71600 -6.00400 23.89600
C  20.84100 -7.42500 25.41400
N  16.84100 -1.78700 25.63700
C  17.26600 -2.56700 24.69400
C  17.04700 -1.94200 23.28700
C  16.12300 -0.73100 23.56000
C  16.07300 -0.74300 25.10600
C  18.39300 -1.56900 22.63000
C  14.74200 -0.86000 22.84700
C  14.78100 -0.68700 21.32000
N  16.09100 -0.30700 28.09200
C  15.35600 0.45700 27.20300
C  14.73300 1.49900 27.97600
C  15.14100 1.35600 29.33600
C  15.94900 0.18700 29.29900
C  13.91400 2.52600 27.44500
C  15.02300 1.80700 30.68600
O  14.43000 2.83800 31.12300
C  15.98100 0.86900 31.51800
C  16.98700 1.71200 32.12800
O  17.99700 2.11700 31.55300
O  16.51500 2.07400 33.36400
C  17.31500 3.14200 34.00000
H  19.94300 -4.54500 30.25900
H  18.26400 -4.34800 24.10300
H  14.78900 0.88500 25.17500
H  18.31000 -0.60500 32.65900
H  18.90200 -3.48800 32.55400
H  20.47800 -0.95900 32.20400
H  21.12500 -2.36600 31.31800
H  21.03200 -2.40200 33.04800
H  16.90400 -3.21300 33.27100
H  15.74800 -2.07700 33.06400
H  16.83100 -2.63200 35.56100
H  16.22600 -1.02300 35.16500
H  21.52400 -6.11200 28.62600
H  20.42400 -7.48300 28.06700
H  20.00400 -6.52900 29.40900
H  20.33800 -8.26100 25.90100
H  21.56700 -6.93500 26.06200
H  21.42700 -7.84100 24.59500
H  16.51200 -2.73300 22.76200
H  16.63000 0.18800 23.26600
H  19.22100 -1.38900 23.31600
H  18.44300 -0.64500 22.05300
H  18.63700 -2.41000 21.98100
H  14.08400 -0.07400 23.21500
H  14.28900 -1.82100 23.09300
H  15.82700 -0.60000 21.02600
H  14.42000 0.31100 21.07100
H  14.22200 -1.38300 20.69500
H  13.25200 2.01200 26.74800
H  14.42200 3.30300 26.87400
H  13.44500 3.05700 28.27400
H  15.32300 0.48100 32.29600
H  17.17000 4.03200 33.38700
H  18.34600 2.83300 34.16900
H  16.87700 3.33200 34.98000


