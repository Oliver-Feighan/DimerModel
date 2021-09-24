%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 52.53200 24.27200 45.27800
C  49.69500 26.14500 43.97100
C  51.04500 21.43900 44.03200
C  55.55800 22.86200 45.08400
C  53.78600 27.15600 46.56500
N  50.65200 23.90000 43.94600
C  49.60000 24.82500 43.60600
C  48.41500 24.05300 43.01900
C  48.70000 22.56700 43.38700
C  50.22600 22.57000 43.84000
C  47.76900 22.00700 44.56400
C  48.04400 24.34300 41.59400
C  49.09500 25.06000 40.71400
H  48.96900 25.39700 39.21700
N  53.20200 22.33000 44.76700
C  52.42300 21.31600 44.25500
C  53.27600 20.09300 44.13100
C  54.57100 20.56600 44.32700
C  54.53600 21.98000 44.69300
C  52.68600 18.78600 43.73500
C  55.85700 19.70500 44.16100
O  55.72900 18.56500 43.75100
C  57.24900 20.21900 44.32400
N  54.40600 25.03400 45.62800
C  55.48900 24.15600 45.60300
C  56.74800 24.82300 46.15700
C  56.16900 26.07400 46.91200
C  54.72900 26.07700 46.39400
C  57.69100 23.98800 47.06100
C  56.92500 27.42300 46.68900
C  57.79700 27.94700 47.92000
N  51.87600 26.25900 45.36900
C  52.48600 27.32200 46.01900
C  51.58600 28.44900 46.02500
C  50.55800 28.02100 45.17400
C  50.79100 26.66600 44.83200
C  51.91200 29.79900 46.68700
C  49.30500 28.37700 44.58200
O  48.71800 29.46600 44.74300
C  48.58400 27.14400 43.87700
C  48.23600 27.64400 42.53200
O  48.98500 28.23500 41.74700
O  46.97500 27.21700 42.27700
C  46.32700 27.76200 41.04300
H  50.57900 20.54200 43.62000
H  56.53800 22.38300 45.02800
H  54.16100 28.03300 47.09600
H  47.65100 24.44400 43.69100
H  48.71700 21.90900 42.51800
H  47.15300 21.26400 44.05900
H  47.21300 22.83200 45.01000
H  48.44700 21.51600 45.26300
H  47.08700 24.86300 41.63300
H  47.84600 23.33500 41.22800
H  49.99900 24.46000 40.82300
H  49.29400 26.05200 41.12000
H  53.06600 18.57900 42.73500
H  51.60100 18.88400 43.75400
H  52.95100 18.03100 44.47600
H  57.60500 20.97300 43.62300
H  58.00300 19.45100 44.15400
H  57.23000 20.50700 45.37500
H  57.21200 25.04800 45.19600
H  56.15200 25.77400 47.96000
H  58.56700 23.60300 46.53800
H  57.16500 23.10200 47.41500
H  57.94800 24.60600 47.92200
H  56.28400 28.25100 46.38500
H  57.54500 27.23200 45.81300
H  58.71300 28.26700 47.42400
H  58.01300 27.19900 48.68300
H  57.20100 28.70700 48.42400
H  51.05000 30.36300 47.04400
H  52.35800 30.33200 45.84800
H  52.47600 29.64100 47.60600
H  47.71400 26.85200 44.46700
H  45.32600 27.34000 40.94700
H  46.88700 27.48700 40.14900
H  46.32300 28.83200 41.24600
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


