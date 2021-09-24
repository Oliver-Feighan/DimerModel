%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg 41.12800 41.85700 26.98600
C  40.12700 44.17900 29.64000
C  41.31000 39.41000 29.48700
C  42.36200 39.82100 24.75900
C  40.70600 44.36300 24.68300
N  40.96700 41.85000 29.29200
C  40.42100 42.87500 30.17300
C  40.33700 42.44400 31.65300
C  40.94600 40.93000 31.52100
C  41.11600 40.67200 30.01100
C  42.22400 40.71400 32.30800
C  38.87900 42.58800 32.18500
C  38.60800 41.92700 33.50300
H  37.80600 42.65600 34.55800
N  41.79000 39.86100 27.08100
C  41.78200 39.01000 28.19200
C  42.26600 37.68700 27.80200
C  42.58400 37.77600 26.38500
C  42.25600 39.19200 25.97700
C  42.31400 36.43000 28.74600
C  43.16500 36.78400 25.44300
O  43.35700 36.97800 24.25000
C  43.44600 35.41200 25.89800
N  41.39000 42.02100 25.04500
C  42.03300 41.09600 24.33000
C  42.18900 41.49200 22.84300
C  41.28100 42.70200 22.76900
C  41.15000 43.12300 24.28400
C  43.62900 41.80100 22.45000
C  39.92300 42.66800 22.11700
C  38.93300 41.74000 22.78800
N  40.64100 43.87900 27.09600
C  40.41100 44.75400 26.01900
C  40.01000 46.03800 26.53900
C  39.83600 45.76900 27.92100
C  40.24100 44.45600 28.22800
C  39.64600 47.28300 25.78400
C  39.34100 46.43200 29.16700
O  38.79000 47.54000 29.34100
C  39.55800 45.41800 30.28800
C  40.41800 46.03500 31.39300
O  41.61900 46.43600 31.36700
O  39.69400 45.95900 32.54600
C  40.38200 46.43100 33.76900
H  41.17000 38.61400 30.22000
H  42.83200 39.17800 24.01200
H  40.57600 45.21400 24.01100
H  40.93200 43.17900 32.19500
H  40.17200 40.21500 31.80000
H  42.06000 39.90200 33.01700
H  42.40500 41.63800 32.85600
H  43.16100 40.69100 31.75200
H  38.23700 42.28000 31.35900
H  38.65800 43.64900 32.30200
H  39.47300 41.65400 34.10700
H  37.99100 41.04400 33.33900
H  41.63300 35.71400 28.28500
H  41.97100 36.63500 29.76000
H  43.34800 36.08500 28.73700
H  43.92500 35.36100 26.87600
H  44.12500 35.02700 25.13800
H  42.45600 34.98600 26.05800
H  41.84700 40.64300 22.25100
H  41.73700 43.52300 22.21500
H  43.81400 41.38100 21.46200
H  44.33400 41.35700 23.15200
H  43.75500 42.87300 22.29400
H  40.06700 42.45700 21.05700
H  39.41200 43.62400 22.23000
H  38.39800 41.12300 22.06700
H  38.14200 42.27600 23.31300
H  39.38000 41.04100 23.49500
H  38.67000 47.63800 26.11400
H  39.60000 47.28100 24.69500
H  40.39300 48.03500 26.03900
H  38.59300 45.21900 30.75500
H  40.75500 45.57700 34.33400
H  39.56300 47.00000 34.20900
H  41.15300 47.17200 33.55400


