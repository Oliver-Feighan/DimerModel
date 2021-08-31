%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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
Mg 46.69300 44.58300 44.38700
C  43.25500 44.72400 43.60300
C  46.82600 41.39700 42.97300
C  50.08000 44.92600 44.22900
C  46.41600 48.03600 45.29800
N  45.21700 43.25300 43.45800
C  43.91500 43.52200 43.22400
C  43.24400 42.37600 42.39900
C  44.33300 41.22100 42.58000
C  45.59600 41.99100 43.08500
C  44.09700 40.05700 43.56900
C  43.05200 42.86300 40.91300
C  41.75200 42.34300 40.23900
H  41.89600 41.90400 38.72500
N  48.30200 43.23100 43.87300
C  48.09600 41.97500 43.36700
C  49.37200 41.38600 43.21500
C  50.37100 42.36600 43.49200
C  49.61800 43.58400 43.94300
C  49.62300 39.97500 42.73800
C  51.83400 42.25100 43.23000
O  52.31500 41.27900 42.67000
C  52.80700 43.26800 43.53300
N  48.03900 46.33200 44.62600
C  49.42400 46.11300 44.66400
C  50.13700 47.33600 45.15000
C  48.99300 48.37700 45.43900
C  47.74400 47.57600 45.11400
C  50.98200 46.96700 46.37600
C  49.23100 49.78300 44.77000
C  49.25400 50.97200 45.68200
N  45.14700 46.13000 44.49400
C  45.21100 47.43000 44.93200
C  43.85400 48.01300 45.01800
C  43.09500 46.96600 44.43100
C  43.89400 45.85700 44.22500
C  43.42700 49.41800 45.41300
C  41.78600 46.66700 43.91400
O  40.77400 47.38200 43.84700
C  41.78000 45.23600 43.37700
C  40.72900 44.42700 44.02200
O  39.86700 43.90600 43.41200
O  40.96200 44.18000 45.34400
C  40.02400 43.15900 45.85900
H  46.89000 40.36900 42.60900
H  51.17100 44.93400 44.24400
H  46.14300 49.07400 45.50200
H  42.28400 42.14500 42.86100
H  44.78800 40.86200 41.65700
H  43.41900 40.36700 44.36500
H  44.98600 39.61600 44.02000
H  43.58300 39.28600 42.99600
H  43.92400 42.61300 40.30800
H  43.01600 43.95000 40.98600
H  41.07000 43.18900 40.32600
H  41.42600 41.56000 40.92400
H  50.58400 39.53700 43.00800
H  49.53400 40.05600 41.65500
H  48.81600 39.25800 42.88900
H  52.77800 43.41800 44.61300
H  52.61000 44.18200 42.97200
H  53.79000 42.90400 43.23400
H  50.70700 47.79400 44.34200
H  48.71000 48.54100 46.47900
H  51.11500 45.89200 46.50000
H  50.52200 47.31200 47.30200
H  52.01900 47.25200 46.20100
H  48.42600 49.87600 44.04100
H  50.21400 49.89500 44.31400
H  48.51600 51.67900 45.30400
H  50.21500 51.48600 45.65500
H  49.14300 50.70900 46.73400
H  43.51500 50.13100 44.59300
H  43.93200 49.78600 46.30700
H  42.35600 49.27300 45.55400
H  41.46600 45.23800 42.33400
H  40.25400 43.00100 46.91300
H  40.12100 42.22900 45.29900
H  38.98500 43.48400 45.90100


