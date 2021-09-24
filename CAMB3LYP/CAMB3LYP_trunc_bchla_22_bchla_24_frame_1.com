%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 9.20100 48.63300 24.80700
C  7.00600 48.65300 27.54500
C  11.64200 49.94100 26.66500
C  11.05000 48.61100 22.15100
C  6.24200 48.21900 22.76700
N  9.22500 49.31100 26.85200
C  8.33600 49.09200 27.82900
C  8.78900 49.67200 29.14400
C  10.10600 50.44700 28.73700
C  10.36100 49.87000 27.35800
C  10.08400 52.01000 28.74700
C  9.14700 48.49500 30.18800
C  8.51500 48.70500 31.60300
H  8.63500 47.36700 32.51500
N  11.09800 48.95000 24.52400
C  12.00800 49.46400 25.43200
C  13.32900 49.31600 24.88600
C  13.24200 48.87600 23.50200
C  11.75200 48.67100 23.33700
C  14.54800 49.78900 25.72900
C  14.23500 48.59100 22.38100
O  13.96800 48.12400 21.25900
C  15.64500 48.82300 22.79500
N  8.62300 48.59400 22.78300
C  9.67100 48.60700 21.89600
C  9.16600 48.63700 20.45300
C  7.69000 48.19700 20.63800
C  7.46000 48.39900 22.13700
C  9.46200 50.04800 19.80600
C  7.33800 46.72100 20.17100
C  5.87200 46.44300 19.80200
N  7.05000 48.42700 25.05400
C  6.04600 48.27100 24.14700
C  4.77900 48.04000 24.80400
C  5.15700 48.11700 26.15700
C  6.52300 48.41500 26.26400
C  3.35900 47.77300 24.25400
C  4.63600 47.99600 27.49000
O  3.59000 47.60400 27.93700
C  5.82500 48.51800 28.45200
C  5.77300 47.51400 29.51800
O  6.28400 46.36600 29.55000
O  5.03000 48.07600 30.48700
C  4.75700 47.29300 31.69300
H  12.47200 50.28300 27.28700
H  11.67900 48.52800 21.26200
H  5.39900 48.00100 22.10800
H  8.05900 50.30600 29.64800
H  11.00500 50.19100 29.29800
H  10.13500 52.34800 27.71200
H  10.98900 52.30100 29.28000
H  9.22300 52.51800 29.18000
H  10.20800 48.27500 30.31200
H  8.82700 47.56600 29.71600
H  7.46000 48.96400 31.51600
H  9.00400 49.45600 32.22300
H  14.89000 48.94400 26.32500
H  14.28100 50.45800 26.54800
H  15.30300 50.38600 25.21800
H  15.98500 48.19800 23.62000
H  15.64000 49.89000 23.02100
H  16.36500 48.70300 21.98500
H  9.75900 47.88800 19.92800
H  7.04700 48.83900 20.03600
H  8.53800 50.61700 19.70400
H  9.90700 50.08500 18.81200
H  10.06100 50.74200 20.39600
H  7.49400 46.14100 21.08100
H  7.92900 46.29300 19.36100
H  5.39900 45.74600 20.49500
H  5.63800 46.05100 18.81200
H  5.34900 47.38100 19.98300
H  3.34300 48.21100 23.25600
H  2.59900 48.22600 24.89200
H  3.25800 46.69000 24.18400
H  5.55100 49.48200 28.88100
H  5.30400 47.73800 32.52400
H  4.92100 46.22700 31.53900
H  3.71000 47.43600 31.96000
Mg -0.18700 43.41100 24.93300
C  1.60700 43.07600 28.07200
C  -2.89800 42.36800 26.64400
C  -1.73200 43.65500 22.14800
C  2.88800 44.15900 23.34500
N  -0.58100 42.80100 27.23000
C  0.26100 42.84200 28.27700
C  -0.50600 42.58400 29.58800
C  -1.80900 41.92300 28.93700
C  -1.80000 42.43200 27.55700
C  -2.07500 40.44900 29.19300
C  -0.70200 43.92400 30.31800
C  -0.86900 43.75600 31.85700
H  -0.56100 42.38900 32.62600
N  -2.11500 43.10200 24.42900
C  -3.10300 42.58200 25.24300
C  -4.30800 42.41600 24.43600
C  -4.00300 42.86400 23.12700
C  -2.55900 43.22400 23.20800
C  -5.57100 41.87400 24.94000
C  -4.95600 43.00700 21.82800
O  -4.58400 43.57000 20.83500
C  -6.46400 42.70700 21.90000
N  0.55900 43.58500 22.97300
C  -0.33400 43.81600 22.01400
C  0.26300 44.22400 20.69400
C  1.62800 44.75300 21.19800
C  1.73900 43.97100 22.55200
C  0.44200 43.02300 19.65300
C  1.76200 46.32800 21.33100
C  0.94100 47.25900 20.39200
N  1.83500 43.48300 25.49800
C  2.95000 43.84600 24.81700
C  4.05200 43.97800 25.65800
C  3.57000 43.63900 26.89000
C  2.22100 43.36000 26.79400
C  5.44000 44.46500 25.25100
C  4.00800 43.52200 28.28900
O  5.10800 43.51800 28.82100
C  2.67500 43.22200 29.10400
C  2.92200 42.05100 29.92700
O  2.92900 40.87200 29.55900
O  3.40200 42.49600 31.10700
C  4.04200 41.44600 31.87900
H  -3.71800 41.86000 27.15600
H  -2.28400 43.68900 21.20700
H  3.79200 44.48700 22.82800
H  -0.02300 41.81500 30.19100
H  -2.65500 42.40200 29.42800
H  -2.95800 40.29300 29.81200
H  -1.21600 39.99100 29.68300
H  -2.18500 39.85300 28.28700
H  -1.65600 44.33500 29.98700
H  0.05400 44.67300 30.08300
H  -1.83500 44.12700 32.20100
H  -0.27300 44.58700 32.23400
H  -5.66200 41.99700 26.01900
H  -5.55700 40.82000 24.66200
H  -6.45800 42.34700 24.51800
H  -7.01900 43.39100 22.54200
H  -6.51800 41.64400 22.13600
H  -6.89500 42.73900 20.89900
H  -0.45600 44.95300 20.32100
H  2.46400 44.44400 20.57100
H  0.05400 43.47000 18.73800
H  -0.13700 42.13000 19.89100
H  1.46800 42.71200 19.46000
H  2.81700 46.59200 21.27300
H  1.31000 46.48200 22.31100
H  -0.04700 47.44600 20.81200
H  0.97500 46.73900 19.43400
H  1.44300 48.20900 20.21200
H  5.25700 45.13800 24.41300
H  6.07800 43.67100 24.86400
H  5.88500 44.97700 26.10400
H  2.48900 44.10200 29.72100
H  4.45500 41.79500 32.82500
H  4.72400 40.78900 31.33900
H  3.30000 40.72800 32.23000


