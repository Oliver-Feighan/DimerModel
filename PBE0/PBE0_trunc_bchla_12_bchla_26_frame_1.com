%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
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
Mg -9.12600 18.73000 43.07000
C  -5.62300 18.44500 42.86700
C  -8.71100 22.08300 42.24600
C  -12.43100 19.03300 42.79300
C  -9.30300 15.26600 42.95700
N  -7.32000 20.02500 42.45300
C  -5.96700 19.76600 42.55400
C  -5.07000 20.98600 42.54900
C  -6.06600 22.15700 42.37300
C  -7.44800 21.42600 42.33500
C  -5.98600 23.25800 43.50700
C  -3.93600 21.03100 41.43000
C  -4.22600 20.35100 40.08200
H  -4.56600 21.25200 38.85800
N  -10.36200 20.30600 42.60500
C  -10.01500 21.60400 42.35900
C  -11.22400 22.42800 42.22000
C  -12.33600 21.54600 42.38200
C  -11.71100 20.20500 42.59800
C  -11.19700 23.93200 41.97500
C  -13.79400 22.10700 42.18800
O  -13.91800 23.28100 41.86800
C  -15.04600 21.14100 42.31600
N  -10.60200 17.36400 42.77700
C  -11.94900 17.71100 42.85200
C  -12.88700 16.37600 42.87100
C  -11.75200 15.22000 42.76300
C  -10.48500 16.00900 42.84500
C  -13.77000 16.28100 44.15000
C  -11.80700 14.33600 41.44400
C  -11.81900 12.80500 41.57400
N  -7.76200 17.11800 43.06200
C  -7.94400 15.76000 43.00400
C  -6.61700 15.11000 42.96500
C  -5.67900 16.17500 43.03100
C  -6.45600 17.36700 43.03800
C  -6.31500 13.64100 42.85700
C  -4.23100 16.45300 42.98800
O  -3.23100 15.72200 43.03800
C  -4.16400 17.97000 43.08300
C  -3.54200 18.42400 44.28100
O  -4.15400 18.45500 45.38000
O  -2.21800 18.84600 44.04800
C  -1.49600 19.49100 45.14000
H  -8.63900 23.16900 42.16100
H  -13.51700 19.11600 42.86600
H  -9.44500 14.18400 42.98300
H  -4.61500 20.94700 43.53800
H  -5.92300 22.50900 41.35200
H  -5.58200 24.18000 43.08900
H  -5.32600 22.85000 44.27300
H  -6.97200 23.48100 43.91600
H  -3.14100 20.46400 41.91600
H  -3.60700 22.06200 41.29800
H  -5.03700 19.62900 40.18300
H  -3.36900 19.72900 39.82700
H  -11.41500 24.23200 40.95000
H  -10.31200 24.38600 42.42100
H  -11.94400 24.43600 42.58800
H  -14.90900 20.56900 43.23400
H  -15.05000 20.49700 41.43700
H  -16.00600 21.63900 42.18000
H  -13.59700 16.50100 42.05400
H  -11.83300 14.61300 43.66500
H  -13.23300 15.66500 44.87100
H  -14.78200 15.93400 43.94100
H  -13.98700 17.24300 44.61400
H  -10.98000 14.66700 40.81500
H  -12.76100 14.64800 41.01900
H  -11.51100 12.31300 40.65100
H  -12.86600 12.60200 41.79800
H  -11.12600 12.60300 42.39100
H  -7.19000 13.22900 42.35500
H  -6.33400 13.32500 43.90000
H  -5.34100 13.39500 42.43300
H  -3.56700 18.21200 42.20400
H  -1.27600 18.88800 46.02100
H  -2.03900 20.37600 45.47300
H  -0.61300 19.94900 44.69400


