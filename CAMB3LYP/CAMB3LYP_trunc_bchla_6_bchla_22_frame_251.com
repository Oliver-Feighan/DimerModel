%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 17.10100 -2.07700 27.77000
C  16.45500 -0.19100 30.58100
C  19.23600 -3.94700 29.68100
C  17.96100 -3.74600 24.96300
C  15.37700 0.23000 25.82200
N  17.83400 -1.91300 29.91400
C  17.31200 -1.24300 30.86800
C  17.86400 -1.55200 32.35400
C  18.99500 -2.58100 31.95700
C  18.76000 -2.80300 30.38800
C  20.50900 -2.04800 32.15900
C  16.79200 -2.13400 33.37400
C  16.96200 -1.78000 34.89400
H  18.25700 -1.14000 35.29400
N  18.31900 -3.63600 27.38700
C  19.10500 -4.38300 28.35000
C  19.74900 -5.47800 27.63200
C  19.48000 -5.36700 26.24500
C  18.51900 -4.24500 26.14900
C  20.50800 -6.46400 28.44500
C  19.91900 -6.24400 25.04900
O  19.71600 -6.00400 23.89600
C  20.84100 -7.42500 25.41400
N  16.84100 -1.78700 25.63700
C  17.26600 -2.56700 24.69400
C  17.04700 -1.94200 23.28700
C  16.12300 -0.73100 23.56000
C  16.07300 -0.74300 25.10600
C  18.39300 -1.56900 22.63000
C  14.74200 -0.86000 22.84700
C  14.78100 -0.68700 21.32000
N  16.09100 -0.30700 28.09200
C  15.35600 0.45700 27.20300
C  14.73300 1.49900 27.97600
C  15.14100 1.35600 29.33600
C  15.94900 0.18700 29.29900
C  13.91400 2.52600 27.44500
C  15.02300 1.80700 30.68600
O  14.43000 2.83800 31.12300
C  15.98100 0.86900 31.51800
C  16.98700 1.71200 32.12800
O  17.99700 2.11700 31.55300
O  16.51500 2.07400 33.36400
C  17.31500 3.14200 34.00000
H  19.94300 -4.54500 30.25900
H  18.26400 -4.34800 24.10300
H  14.78900 0.88500 25.17500
H  18.31000 -0.60500 32.65900
H  18.90200 -3.48800 32.55400
H  20.47800 -0.95900 32.20400
H  21.12500 -2.36600 31.31800
H  21.03200 -2.40200 33.04800
H  16.90400 -3.21300 33.27100
H  15.74800 -2.07700 33.06400
H  16.83100 -2.63200 35.56100
H  16.22600 -1.02300 35.16500
H  21.52400 -6.11200 28.62600
H  20.42400 -7.48300 28.06700
H  20.00400 -6.52900 29.40900
H  20.33800 -8.26100 25.90100
H  21.56700 -6.93500 26.06200
H  21.42700 -7.84100 24.59500
H  16.51200 -2.73300 22.76200
H  16.63000 0.18800 23.26600
H  19.22100 -1.38900 23.31600
H  18.44300 -0.64500 22.05300
H  18.63700 -2.41000 21.98100
H  14.08400 -0.07400 23.21500
H  14.28900 -1.82100 23.09300
H  15.82700 -0.60000 21.02600
H  14.42000 0.31100 21.07100
H  14.22200 -1.38300 20.69500
H  13.25200 2.01200 26.74800
H  14.42200 3.30300 26.87400
H  13.44500 3.05700 28.27400
H  15.32300 0.48100 32.29600
H  17.17000 4.03200 33.38700
H  18.34600 2.83300 34.16900
H  16.87700 3.33200 34.98000
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


