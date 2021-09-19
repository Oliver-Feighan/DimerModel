%nproc=24
%mem=175gb
#p wB97XD/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
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
Mg 0.26000 43.34200 24.75600
C  2.25300 43.09900 27.61800
C  -2.44200 42.42000 26.69400
C  -1.57400 42.87700 21.91500
C  3.14600 43.60500 22.86000
N  -0.01100 42.92700 26.88700
C  0.88400 43.10000 27.86800
C  0.29200 43.09400 29.28500
C  -1.20300 42.68500 28.99500
C  -1.29800 42.75900 27.37000
C  -1.68200 41.32300 29.54200
C  0.44600 44.41100 30.04600
C  0.46000 44.23200 31.57500
H  0.30000 42.82500 32.16800
N  -1.74500 42.84700 24.34000
C  -2.74500 42.55800 25.29700
C  -3.99500 42.34000 24.59700
C  -3.74800 42.40500 23.17500
C  -2.33100 42.75300 23.08200
C  -5.28700 42.18000 25.35000
C  -4.75800 42.32200 21.94400
O  -4.44000 42.53400 20.79900
C  -6.18400 41.93600 22.27200
N  0.72800 43.16000 22.71000
C  -0.18900 43.01100 21.71100
C  0.47200 43.04900 20.28400
C  1.92800 43.34600 20.63900
C  1.92700 43.36400 22.15600
C  0.26000 41.87700 19.24300
C  2.38800 44.69400 20.03600
C  2.04300 46.00200 20.78200
N  2.32100 43.50900 25.09100
C  3.38200 43.61100 24.22900
C  4.61600 43.69300 24.93600
C  4.22700 43.48400 26.33500
C  2.89200 43.36300 26.33900
C  5.99200 43.88800 24.47400
C  4.71200 43.25500 27.68900
O  5.85600 43.34100 28.14800
C  3.47500 43.00100 28.60000
C  3.56600 41.71300 29.28000
O  3.33400 40.61600 28.78300
O  4.17100 41.89000 30.49900
C  4.75700 40.73100 31.20000
H  -3.23100 42.10200 27.37900
H  -2.18100 42.91300 21.00800
H  4.04900 43.78500 22.27300
H  0.67700 42.26400 29.87700
H  -1.92300 43.45600 29.26700
H  -2.57300 41.63600 30.08600
H  -0.94100 40.89700 30.21800
H  -1.87200 40.62800 28.72400
H  -0.32500 45.13200 29.77700
H  1.38300 44.89600 29.77100
H  -0.37400 44.87500 31.85700
H  1.32100 44.72000 32.03400
H  -5.00100 42.06900 26.39600
H  -5.81000 41.30000 24.97500
H  -5.96200 43.02800 25.22500
H  -6.37400 40.88800 22.50600
H  -6.81500 42.27400 21.45000
H  -6.38700 42.55900 23.14300
H  -0.03500 43.91100 19.85100
H  2.51500 42.52100 20.23700
H  1.17400 41.51200 18.77600
H  -0.34600 42.26100 18.42300
H  -0.30900 41.09200 19.74200
H  1.97300 44.71300 19.02900
H  3.46900 44.59400 19.93400
H  1.75200 46.69500 19.99300
H  2.89700 46.44600 21.29300
H  1.29400 45.88900 21.56500
H  6.47500 44.58100 25.16200
H  6.09900 44.26000 23.45500
H  6.57000 42.97800 24.63600
H  3.48200 43.84300 29.29200
H  4.58400 39.79600 30.66600
H  4.26600 40.64900 32.17000
H  5.80700 41.02000 31.22900


