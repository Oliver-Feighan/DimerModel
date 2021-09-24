%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 16.95800 -2.04300 28.05800
C  15.67900 0.03700 30.54100
C  18.78700 -3.72100 30.44000
C  18.31300 -3.78900 25.50500
C  15.32700 0.02400 25.64200
N  17.25900 -1.76300 30.19400
C  16.62800 -0.88400 31.01900
C  16.98400 -1.14100 32.48900
C  18.10000 -2.26900 32.35300
C  18.02400 -2.61400 30.86200
C  19.48200 -1.75800 32.74200
C  15.83300 -1.72700 33.38200
C  15.87500 -1.41800 34.85600
H  17.10800 -0.90300 35.45200
N  18.26300 -3.57700 28.00200
C  18.92400 -4.12800 29.11800
C  19.62000 -5.29300 28.63800
C  19.47300 -5.38800 27.23100
C  18.65200 -4.21500 26.84200
C  20.55500 -6.11300 29.50100
C  19.94500 -6.39800 26.27100
O  19.91000 -6.24200 25.06900
C  20.68300 -7.69100 26.77500
N  16.72600 -2.04100 25.82800
C  17.40600 -2.83200 25.05300
C  17.12300 -2.53800 23.54900
C  16.12900 -1.39100 23.56700
C  16.05500 -1.03300 25.09500
C  18.49500 -2.31000 22.75600
C  14.77700 -1.75000 22.91500
C  14.57500 -1.34500 21.40700
N  15.59300 -0.54700 28.03400
C  15.01200 0.20800 26.98500
C  14.13500 1.19500 27.48900
C  14.45500 1.19400 28.91300
C  15.31700 0.08300 29.13900
C  13.41400 2.25000 26.67000
C  14.16900 1.92500 30.09200
O  13.43900 2.88400 30.30900
C  14.95100 1.07800 31.29700
C  15.78000 2.10000 32.01800
O  16.65500 2.79200 31.48700
O  15.36200 2.29400 33.36800
C  16.13600 3.30900 34.14900
H  19.34000 -4.20600 31.24800
H  18.93000 -4.26900 24.74300
H  14.76600 0.65500 24.95000
H  17.39400 -0.20500 32.86900
H  17.90900 -3.13400 32.98800
H  19.41600 -0.81500 33.28400
H  20.16500 -1.70100 31.89400
H  19.97000 -2.35300 33.51400
H  15.82100 -2.80600 33.23200
H  14.89600 -1.40800 32.92500
H  15.72700 -2.36200 35.38000
H  15.08500 -0.70500 35.09200
H  21.59900 -6.03400 29.19800
H  20.20800 -7.14700 29.48400
H  20.60300 -5.83600 30.55400
H  21.71600 -7.39100 26.95300
H  20.76400 -8.41600 25.96500
H  20.25400 -8.04400 27.71200
H  16.67800 -3.47000 23.20000
H  16.55100 -0.48300 23.13900
H  18.49400 -1.34000 22.25900
H  18.63800 -3.00100 21.92500
H  19.26300 -2.34000 23.53000
H  14.01500 -1.24300 23.50800
H  14.64200 -2.82300 23.05200
H  14.22000 -0.31400 21.38800
H  13.95900 -1.97800 20.76800
H  15.61200 -1.30800 21.07500
H  12.51000 1.73500 26.34600
H  13.99600 2.54300 25.79700
H  13.31000 3.23400 27.12900
H  14.20300 0.59900 31.93000
H  16.03100 4.30700 33.72500
H  17.20500 3.14900 34.29200
H  15.57400 3.26700 35.08200
Mg -9.66500 18.33400 42.78400
C  -6.24000 17.80500 42.49600
C  -9.09900 21.67600 42.05300
C  -12.91200 18.87400 42.48300
C  -10.14800 14.91400 43.00500
N  -7.91600 19.56100 42.17100
C  -6.59800 19.16200 42.25800
C  -5.63300 20.29600 42.07600
C  -6.57300 21.55800 42.31900
C  -7.95100 20.93600 42.13600
C  -6.37100 22.54100 43.56800
C  -4.81400 20.28500 40.70300
C  -5.51800 19.69700 39.44400
H  -5.27900 20.48100 38.17600
N  -10.85700 20.01500 42.21400
C  -10.42000 21.32200 41.96500
C  -11.51000 22.18000 41.58900
C  -12.66600 21.36900 41.68400
C  -12.20600 20.01400 42.11100
C  -11.29500 23.68900 41.46200
C  -14.07900 21.93000 41.40500
O  -14.25800 23.08100 40.99400
C  -15.20100 21.03000 41.71200
N  -11.24800 17.01300 42.63400
C  -12.52300 17.52900 42.79400
C  -13.59000 16.43700 43.02800
C  -12.64200 15.17100 42.94600
C  -11.28500 15.70200 42.80600
C  -14.41600 16.47100 44.27500
C  -12.97300 14.21900 41.68800
C  -13.02000 12.77100 42.10100
N  -8.45700 16.64500 42.95800
C  -8.80300 15.30900 43.04100
C  -7.58100 14.47400 43.08900
C  -6.53000 15.43300 42.87800
C  -7.12900 16.69600 42.80500
C  -7.42000 13.00900 43.03700
C  -5.08700 15.65900 42.64900
O  -4.19400 14.85200 42.53200
C  -4.84400 17.18800 42.43700
C  -3.86600 17.68600 43.40500
O  -2.75700 18.18900 43.22300
O  -4.28400 17.22800 44.68200
C  -3.45300 17.67600 45.82000
H  -8.86100 22.74200 42.07500
H  -13.99400 18.96100 42.60100
H  -10.40600 13.86400 43.15400
H  -4.85200 20.25200 42.83500
H  -6.35600 22.13700 41.42100
H  -5.49500 22.34800 44.18700
H  -7.18200 22.50000 44.29500
H  -6.14100 23.54200 43.20200
H  -3.84000 19.81700 40.84200
H  -4.55300 21.32400 40.50200
H  -6.59000 19.59400 39.61500
H  -5.14100 18.69400 39.24600
H  -11.94900 24.16700 40.73200
H  -10.29700 23.87500 41.06600
H  -11.49900 24.10000 42.45100
H  -15.27200 20.92300 42.79400
H  -15.05600 20.07500 41.20600
H  -16.08400 21.54300 41.33100
H  -14.18800 16.36800 42.11900
H  -12.69300 14.77700 43.96100
H  -14.92900 15.51000 44.31600
H  -15.20300 17.20800 44.11500
H  -13.79200 16.69500 45.14100
H  -12.21000 14.42700 40.93800
H  -13.86700 14.45400 41.10900
H  -12.15800 12.27100 41.66000
H  -13.81800 12.22200 41.60000
H  -12.89800 12.69100 43.18100
H  -8.34400 12.43000 43.00100
H  -6.93100 12.73100 43.97000
H  -6.96900 12.89900 42.05100
H  -4.46800 17.26400 41.41700
H  -4.21000 17.98200 46.54300
H  -2.74000 18.49300 45.70900
H  -3.06600 16.75200 46.25000


