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
Mg 47.67100 15.89000 28.45400
C  45.34400 15.65500 31.04800
C  49.48000 18.07500 30.51500
C  49.68700 16.65900 25.89600
C  45.30600 14.62100 26.30000
N  47.50700 16.71700 30.64000
C  46.40900 16.47700 31.40500
C  46.52000 17.23500 32.78700
C  47.75800 18.22000 32.48000
C  48.27100 17.69900 31.12700
C  47.16600 19.65600 32.25100
C  46.90400 16.16800 33.87700
C  47.20400 16.68600 35.31000
H  48.31600 16.29900 36.23900
N  49.43000 17.01400 28.32000
C  50.06900 17.82700 29.27100
C  51.20000 18.39100 28.66300
C  51.38300 17.85100 27.39700
C  50.14500 17.12800 27.15700
C  52.11500 19.27100 29.45000
C  52.51500 18.07600 26.34100
O  52.57400 17.52300 25.24800
C  53.78800 18.86100 26.66900
N  47.44200 15.76700 26.39600
C  48.48400 16.01700 25.57200
C  48.26000 15.53000 24.12700
C  46.80100 15.00400 24.22400
C  46.51000 15.00800 25.71400
C  48.61100 16.58300 23.02400
C  46.48600 13.65000 23.53100
C  46.96600 12.43800 24.21000
N  45.75600 15.16700 28.59900
C  44.93800 14.69700 27.68600
C  43.69800 14.18000 28.27400
C  43.86300 14.56100 29.64400
C  45.08900 15.18100 29.78600
C  42.77800 13.18000 27.74800
C  43.21100 14.49100 30.94000
O  42.19200 13.93500 31.28200
C  44.14500 15.27600 31.86900
C  43.41900 16.43800 32.55500
O  43.17000 17.49000 32.02600
O  43.04200 16.05500 33.82000
C  42.33900 17.13300 34.55300
H  50.09800 18.64700 31.20900
H  50.41400 16.70200 25.08200
H  44.81200 14.14000 25.45400
H  45.63400 17.83000 33.01000
H  48.52200 18.24600 33.25800
H  47.13800 20.01300 31.22200
H  47.84000 20.42900 32.62100
H  46.16600 19.79300 32.66200
H  47.75800 15.61000 33.49400
H  46.09800 15.43400 33.89800
H  46.24800 16.58700 35.82500
H  47.42100 17.73200 35.09700
H  51.80300 19.48300 30.47300
H  52.34700 20.14300 28.83700
H  53.02500 18.67700 29.53900
H  53.61800 19.92500 26.83400
H  54.61200 18.70300 25.97300
H  54.19200 18.41700 27.57900
H  49.00300 14.73300 24.13500
H  46.20400 15.80500 23.78800
H  49.40000 17.25400 23.36200
H  47.72600 17.18100 22.80500
H  48.98600 16.12600 22.10800
H  46.99700 13.77600 22.57700
H  45.42300 13.59700 23.29500
H  46.25100 12.14600 24.97900
H  47.97200 12.58100 24.60400
H  46.98800 11.57800 23.54000
H  43.37600 12.59600 27.04800
H  41.93300 13.49000 27.13300
H  42.41600 12.61300 28.60500
H  44.52600 14.58100 32.61800
H  41.57800 16.67700 35.18600
H  41.70900 17.79000 33.95300
H  43.02300 17.66800 35.21100


