%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 40.61000 41.38800 26.95600
C  39.86800 43.56800 29.58900
C  41.56600 39.09100 29.09600
C  42.01100 39.61800 24.36600
C  40.07100 44.02500 24.67100
N  40.59100 41.34900 29.06300
C  40.20000 42.29100 29.98900
C  40.16900 41.71800 31.40000
C  40.93300 40.36900 31.20600
C  40.94200 40.20500 29.72600
C  42.27000 40.16100 31.88400
C  38.74600 41.61400 32.03600
C  38.70900 41.06400 33.46600
H  37.95300 41.90400 34.48500
N  41.75400 39.53700 26.76500
C  41.90400 38.67000 27.78800
C  42.49400 37.45600 27.29600
C  42.73100 37.64100 25.92300
C  42.16600 38.98800 25.63800
C  42.83400 36.26800 28.17100
C  43.30000 36.69700 24.91500
O  43.61700 36.97700 23.72800
C  43.50600 35.23000 25.27200
N  41.00200 41.78000 24.87100
C  41.64700 40.91600 24.05300
C  41.55500 41.44600 22.52200
C  40.51800 42.61700 22.64600
C  40.44500 42.80800 24.12700
C  42.99600 41.87100 22.03700
C  39.09100 42.15200 22.17700
C  38.49200 40.88200 22.76100
N  40.00500 43.44200 27.05400
C  39.84200 44.33300 26.02900
C  39.43100 45.53300 26.59500
C  39.42300 45.36800 28.01000
C  39.78400 44.06700 28.24300
C  39.15600 46.83600 25.93100
C  39.20700 45.98800 29.33700
O  38.93200 47.12500 29.70300
C  39.39600 44.78700 30.41700
C  40.24300 45.23600 31.55400
O  41.22600 46.02900 31.45500
O  39.66700 44.79200 32.66300
C  40.07000 45.50000 33.86300
H  41.79000 38.25100 29.75700
H  42.42000 39.05300 23.52500
H  39.79000 44.77400 23.92900
H  40.67200 42.47500 32.00200
H  40.36400 39.50500 31.55100
H  42.49900 40.80900 32.73000
H  43.06200 40.23200 31.13800
H  42.35700 39.09400 32.09000
H  38.09900 41.02000 31.39000
H  38.28900 42.60300 32.04400
H  39.66400 41.20800 33.97100
H  38.39800 40.01900 33.48100
H  42.71500 36.43800 29.24100
H  43.88400 35.97700 28.18900
H  42.26200 35.37200 27.92800
H  44.31400 35.17500 26.00100
H  43.97700 34.67500 24.46100
H  42.58300 34.72100 25.55200
H  41.23200 40.63200 21.87300
H  40.99300 43.46400 22.15100
H  43.79600 41.51800 22.68800
H  43.14700 42.94900 21.97800
H  43.27300 41.44900 21.07200
H  39.29200 42.04500 21.11100
H  38.30900 42.90300 22.28500
H  39.25300 40.32600 23.30900
H  37.94200 40.24300 22.07100
H  37.70100 41.27300 23.40100
H  38.42700 46.78300 25.12200
H  40.10600 47.12900 25.48500
H  38.86200 47.61200 26.63700
H  38.38700 44.62200 30.79400
H  39.19900 45.70100 34.48800
H  40.53400 46.47500 33.71600
H  40.74300 44.82600 34.39400
Mg 8.73000 48.00700 24.72600
C  6.54700 48.35400 27.54300
C  11.38500 48.83200 26.74200
C  10.78900 48.16000 22.05300
C  5.85600 48.14800 22.68600
N  8.91100 48.51400 26.88400
C  7.92400 48.53400 27.88200
C  8.57600 48.79800 29.21500
C  9.96600 49.43600 28.79200
C  10.09000 48.94100 27.32100
C  9.83300 50.99700 28.89900
C  8.62100 47.55100 30.16700
C  8.25400 47.92900 31.51300
H  7.91700 46.69900 32.49700
N  10.84000 48.27700 24.55100
C  11.74400 48.56500 25.49500
C  13.05300 48.54500 24.91800
C  12.92300 48.37200 23.55600
C  11.42000 48.19500 23.29600
C  14.33000 48.90700 25.74100
C  14.02400 48.24600 22.48300
O  13.70400 48.13100 21.32700
C  15.48100 48.17400 22.87500
N  8.33400 48.41200 22.66100
C  9.38600 48.28200 21.76100
C  8.87500 48.15200 20.36300
C  7.28100 48.12500 20.55900
C  7.12600 48.16700 22.06100
C  9.40100 49.20400 19.44600
C  6.46300 46.89200 19.88400
C  5.44500 47.20900 18.86100
N  6.61400 48.15800 24.98900
C  5.61400 48.20500 24.04400
C  4.27700 48.23000 24.78200
C  4.59000 48.17400 26.15500
C  6.02700 48.18500 26.22900
C  2.90700 48.32700 24.09700
C  3.98700 48.26800 27.45000
O  2.78100 48.32100 27.78100
C  5.21700 48.56000 28.40400
C  5.00500 47.63900 29.64300
O  5.64000 46.57500 29.76100
O  4.15000 48.22900 30.51400
C  3.80400 47.24100 31.59300
H  12.25300 48.87100 27.40400
H  11.43600 48.18800 21.17300
H  4.90900 48.09600 22.14500
H  8.00400 49.60800 29.66700
H  10.77600 49.10700 29.44200
H  9.23900 51.44000 28.10000
H  10.83100 51.43100 28.94300
H  9.35800 51.24300 29.84900
H  9.67300 47.28200 30.26500
H  8.05700 46.70900 29.76500
H  7.29600 48.44900 31.52800
H  9.07600 48.51700 31.92200
H  14.96900 48.03700 25.59000
H  14.16900 49.02000 26.81300
H  14.74800 49.84700 25.38100
H  16.11600 48.17200 21.98900
H  15.42200 47.21200 23.38400
H  15.69700 49.01800 23.53000
H  9.23500 47.16500 20.07200
H  6.89900 49.04100 20.10900
H  10.07700 48.66800 18.78100
H  10.03800 49.93000 19.95100
H  8.64800 49.76300 18.89000
H  5.91200 46.49600 20.73700
H  7.19200 46.14700 19.56800
H  4.44900 46.91600 19.19200
H  5.77800 46.47600 18.12700
H  5.55100 48.21100 18.44300
H  2.66500 49.38900 24.07600
H  2.13000 47.78700 24.63800
H  2.97200 47.86900 23.11000
H  5.09400 49.61500 28.64500
H  3.12300 46.44700 31.28700
H  3.39400 47.68200 32.50100
H  4.75700 46.82700 31.92300


