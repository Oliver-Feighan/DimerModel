%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 2.79300 -0.99900 44.50000
C  6.05500 -0.17400 43.38400
C  1.60600 1.42400 42.43900
C  -0.30000 -2.33900 44.94400
C  4.15400 -3.74300 46.12200
N  3.83100 0.31200 42.98000
C  5.15200 0.58300 42.74000
C  5.32600 1.71500 41.79700
C  3.88200 2.40600 41.78900
C  3.05100 1.33800 42.41400
C  3.87900 3.72200 42.55000
C  5.83700 1.13000 40.47100
C  7.17200 1.73800 39.86600
H  7.22000 2.30200 38.41400
N  0.88600 -0.54700 43.86000
C  0.60800 0.57500 43.06700
C  -0.81500 0.72700 43.01300
C  -1.42800 -0.37100 43.70600
C  -0.26400 -1.10400 44.27200
C  -1.44600 1.92200 42.31000
C  -2.89700 -0.69100 43.87000
O  -3.66200 0.20700 43.41300
C  -3.55500 -1.99100 44.41300
N  1.99700 -2.80800 45.42300
C  0.68200 -3.10900 45.50100
C  0.40600 -4.46500 46.20300
C  1.84700 -4.85700 46.68000
C  2.74500 -3.73600 45.97300
C  -0.59100 -4.26300 47.33000
C  2.30300 -6.22900 46.28700
C  2.45600 -7.17500 47.44700
N  4.64000 -1.73900 44.89100
C  5.11000 -2.84400 45.61200
C  6.47200 -2.80600 45.78400
C  6.89300 -1.74700 44.95400
C  5.75400 -1.18800 44.37900
C  7.21700 -3.66700 46.78200
C  8.08600 -0.96000 44.54300
O  9.22400 -1.08600 44.87400
C  7.50800 0.01200 43.44600
C  8.17400 1.35500 43.65200
O  9.07700 1.80400 42.94600
O  7.57500 2.04000 44.73700
C  8.09900 3.39600 44.89800
H  1.22700 2.33000 41.96100
H  -1.28000 -2.75100 45.19600
H  4.63400 -4.67900 46.41300
H  6.07700 2.40100 42.18900
H  3.61200 2.61400 40.75400
H  4.50800 3.70600 43.43900
H  2.87200 3.89300 42.93200
H  4.11600 4.58800 41.93200
H  5.03800 1.25300 39.74000
H  5.94100 0.04500 40.47900
H  7.81400 0.85800 39.91300
H  7.65200 2.40300 40.58300
H  -1.04200 2.04400 41.30500
H  -1.16400 2.79800 42.89400
H  -2.51800 1.88200 42.11200
H  -3.31300 -1.96500 45.47500
H  -3.08800 -2.84600 43.92400
H  -4.62700 -2.03100 44.21800
H  -0.01900 -5.22600 45.54900
H  1.95500 -4.64800 47.74400
H  -1.44700 -3.72600 46.92100
H  -0.11100 -3.66200 48.10300
H  -0.92100 -5.12200 47.91300
H  3.28200 -6.24400 45.80800
H  1.69000 -6.76000 45.55900
H  1.86400 -8.06400 47.23300
H  2.03500 -6.78100 48.37200
H  3.49000 -7.42400 47.68800
H  8.18600 -3.18800 46.92600
H  7.36400 -4.68500 46.42200
H  6.59100 -3.65900 47.67400
H  7.98500 -0.36600 42.54100
H  7.35600 4.03900 45.37000
H  8.23900 3.88500 43.93500
H  9.04800 3.42700 45.43400
Mg 47.33900 34.91800 28.07800
C  45.68800 33.19800 30.68300
C  48.10100 37.44900 30.24400
C  48.24300 36.72400 25.47200
C  46.07200 32.47700 25.83700
N  47.02000 35.19800 30.25900
C  46.36000 34.38500 31.06000
C  46.67100 34.77400 32.49400
C  47.03300 36.30500 32.38100
C  47.41200 36.35800 30.85400
C  45.84200 37.21500 32.86800
C  47.88800 34.00400 33.06800
C  47.98300 33.94800 34.65600
H  46.95600 34.64300 35.55400
N  48.24300 36.70300 27.90700
C  48.56900 37.63000 28.93000
C  49.30300 38.76000 28.31400
C  49.26700 38.57700 26.88500
C  48.57000 37.29400 26.64900
C  49.93800 39.85800 29.09200
C  49.83100 39.42600 25.68000
O  49.78800 39.06500 24.49600
C  50.32700 40.81500 25.96000
N  47.07800 34.69300 25.91700
C  47.58500 35.54400 25.06500
C  47.43000 35.17200 23.61700
C  46.71000 33.78100 23.71200
C  46.68600 33.56800 25.25900
C  46.76800 36.29200 22.80900
C  47.43400 32.63900 23.01000
C  46.53900 31.60900 22.23200
N  46.16500 33.09500 28.16200
C  45.77400 32.19200 27.15400
C  45.18200 31.05300 27.76700
C  44.95200 31.44600 29.12400
C  45.59700 32.64100 29.31900
C  44.92100 29.73400 27.09800
C  44.48300 31.08700 30.48400
O  43.89600 30.09800 30.87100
C  44.73800 32.34900 31.43800
C  44.99200 31.90700 32.84200
O  45.89000 31.21500 33.22200
O  43.96300 32.32800 33.70500
C  43.67300 31.45600 34.86300
H  48.20800 38.37100 30.82000
H  48.63100 37.17100 24.55500
H  45.76300 31.70600 25.12800
H  45.78600 34.43600 33.03300
H  47.91000 36.55500 32.97900
H  44.99600 36.59900 33.17100
H  45.50600 37.70100 31.95200
H  46.07700 37.86100 33.71400
H  48.75200 34.55000 32.68700
H  47.86600 33.00400 32.63500
H  48.97800 34.33500 34.88100
H  48.01300 32.88400 34.89200
H  50.02200 39.76300 30.17400
H  49.40400 40.80200 28.98700
H  50.95300 39.89900 28.69600
H  51.37600 40.86200 26.25200
H  49.62900 41.09200 26.75100
H  50.11900 41.39300 25.06000
H  48.46800 35.04700 23.30900
H  45.69100 33.78300 23.32700
H  46.13600 36.95200 23.40400
H  46.11300 35.78300 22.10300
H  47.49100 36.91500 22.28200
H  48.02000 32.05200 23.71700
H  48.21400 33.13900 22.43600
H  45.47900 31.78000 22.41900
H  46.81500 30.64200 22.65200
H  46.65400 31.52300 21.15100
H  44.41500 29.99000 26.16700
H  44.28000 29.05700 27.66300
H  45.86600 29.26100 26.82900
H  43.79600 32.88200 31.56300
H  43.90900 32.01600 35.76800
H  44.22800 30.51800 34.82500
H  42.60300 31.25600 34.80700


