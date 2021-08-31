%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 8.82800 2.97900 28.39100
C  10.36500 1.60400 31.10700
C  7.59700 5.54500 30.49200
C  7.11400 4.23300 25.80200
C  10.13400 0.47400 26.43700
N  8.69300 3.35100 30.57800
C  9.52500 2.68500 31.48700
C  9.40600 3.27200 32.92200
C  8.46000 4.46400 32.61800
C  8.16100 4.43600 31.13000
C  7.12300 4.69600 33.44300
C  10.76000 3.72700 33.62500
C  10.50900 4.32800 35.03400
H  10.96800 3.41400 36.17100
N  7.54400 4.72800 28.17300
C  7.38800 5.73600 29.10000
C  6.47100 6.69300 28.49200
C  6.33200 6.29300 27.12400
C  7.02200 5.05300 26.95500
C  5.95900 7.85000 29.26000
C  5.41800 7.05100 26.09900
O  5.12000 6.60100 24.97400
C  4.78500 8.36300 26.48300
N  8.63500 2.44800 26.51900
C  7.83100 3.12600 25.56000
C  7.67200 2.29800 24.30500
C  8.94400 1.45600 24.38400
C  9.32300 1.47600 25.85800
C  6.36700 1.43300 24.30400
C  10.11600 2.06600 23.58300
C  10.22200 1.75100 22.11500
N  9.99100 1.29200 28.61900
C  10.57300 0.45400 27.74600
C  11.35400 -0.51100 28.46200
C  11.33200 -0.10700 29.78000
C  10.51600 1.03400 29.82700
C  12.08500 -1.61500 27.78100
C  11.82500 -0.42000 31.05300
O  12.58200 -1.29100 31.43600
C  11.20500 0.70900 32.00300
C  10.44000 0.02800 33.08200
O  9.42100 -0.60800 32.96400
O  10.93800 0.32100 34.28700
C  10.21200 -0.27100 35.47600
H  7.38800 6.37900 31.16500
H  6.56600 4.63400 24.94700
H  10.68600 -0.27500 25.86600
H  8.89500 2.54300 33.55100
H  9.00500 5.40000 32.74400
H  6.85500 5.73300 33.64700
H  7.10100 4.02900 34.30500
H  6.33500 4.37900 32.76000
H  11.20100 4.41000 32.90000
H  11.41200 2.85700 33.70200
H  9.47300 4.65600 35.11900
H  11.13500 5.21700 35.10500
H  4.89200 7.69700 29.42300
H  6.27900 8.73500 28.71000
H  6.46700 7.93200 30.22100
H  4.23800 8.77300 25.63400
H  5.63600 9.03200 26.61800
H  4.16200 8.37000 27.37700
H  7.64200 2.93000 23.41700
H  8.80000 0.42300 24.06800
H  6.45900 0.58300 24.97900
H  6.13800 1.09500 23.29300
H  5.57600 2.10700 24.63100
H  11.08800 1.85700 24.02800
H  10.10000 3.15100 23.68700
H  10.17300 2.70400 21.58900
H  9.27100 1.26300 21.90100
H  11.07000 1.12100 21.84600
H  13.12700 -1.33900 27.62400
H  11.63900 -2.06200 26.89200
H  12.10800 -2.42700 28.50800
H  12.02700 1.28100 32.43400
H  9.15600 -0.02800 35.35600
H  10.52300 0.16000 36.42800
H  10.37400 -1.34900 35.46000
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


