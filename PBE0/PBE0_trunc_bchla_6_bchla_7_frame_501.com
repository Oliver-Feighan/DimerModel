%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg 25.26700 0.12000 29.48500
C  27.14100 -0.08700 32.36300
C  22.44200 0.65100 31.48400
C  23.33900 0.24400 26.71700
C  27.97000 -0.76800 27.63500
N  24.82700 -0.00300 31.72000
C  25.83000 0.13900 32.68700
C  25.26100 0.57400 34.05200
C  23.73000 0.68200 33.70700
C  23.61100 0.48500 32.21100
C  22.87000 -0.30200 34.49000
C  25.80700 1.87500 34.74000
C  26.64900 1.68700 36.11600
H  25.98700 2.43400 37.27800
N  23.21300 0.29800 29.14100
C  22.23900 0.64000 30.09300
C  20.97400 1.01400 29.45200
C  21.16000 0.82900 27.97700
C  22.63600 0.50000 27.85200
C  19.74900 1.44500 30.22800
C  20.23300 0.96600 26.86300
O  20.59800 0.66100 25.71800
C  18.79000 1.40700 27.03500
N  25.61300 -0.29200 27.55300
C  24.69700 -0.04800 26.51300
C  25.29200 -0.18000 25.05300
C  26.81800 -0.54700 25.43400
C  26.79800 -0.51500 26.97900
C  24.57500 -1.29000 24.13700
C  27.83300 0.46900 24.80600
C  29.13600 -0.06400 24.29800
N  27.18000 -0.44400 29.94400
C  28.21700 -0.73500 29.05300
C  29.45900 -0.85000 29.73800
C  29.08500 -0.49600 31.08700
C  27.69300 -0.34600 31.12700
C  30.82000 -1.11200 29.24500
C  29.58400 -0.41400 32.42400
O  30.70000 -0.50100 32.94800
C  28.28000 -0.28500 33.33500
C  28.18400 -1.48700 34.20500
O  27.52600 -2.44100 34.06000
O  28.92900 -1.30000 35.36700
C  28.48600 -2.07100 36.51000
H  21.58300 0.90400 32.10900
H  22.84400 0.32400 25.74700
H  28.87000 -0.90500 27.03200
H  25.38700 -0.24700 34.75700
H  23.31500 1.64400 34.00900
H  22.43700 -1.01900 33.79300
H  22.09100 0.33900 34.90100
H  23.35600 -0.79100 35.33400
H  24.95600 2.52500 34.94500
H  26.49900 2.40700 34.08700
H  27.58200 2.20000 35.88200
H  26.76100 0.62900 36.35200
H  19.34600 0.54600 30.69600
H  19.02800 1.84900 29.51900
H  19.97500 2.10600 31.06500
H  18.28300 1.10300 26.11900
H  18.67500 2.47500 27.22300
H  18.18400 0.86400 27.76000
H  25.32400 0.80000 24.57800
H  27.16700 -1.50900 25.05800
H  25.34100 -1.79600 23.54900
H  23.91500 -0.92500 23.35100
H  24.00900 -2.02000 24.71600
H  28.05800 1.32200 25.44600
H  27.45300 0.92700 23.89300
H  29.88300 0.35500 24.97200
H  29.54600 0.28700 23.35100
H  29.13600 -1.15300 24.35500
H  30.88400 -0.76300 28.21400
H  30.89900 -2.19900 29.25700
H  31.53400 -0.73600 29.97800
H  28.35300 0.59900 33.96800
H  28.85600 -3.07500 36.30500
H  27.39700 -2.05300 36.47000
H  28.92200 -1.70000 37.43800


