%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 40.54200 8.71200 29.57000
C  42.46900 10.35300 31.83700
C  38.24000 7.70200 31.93300
C  38.87500 6.96300 27.14400
C  43.10200 9.33000 27.08900
N  40.50700 8.83000 31.66000
C  41.28100 9.79100 32.39300
C  40.70700 10.01900 33.73400
C  39.62800 8.84800 33.81600
C  39.41100 8.39300 32.39000
C  40.09700 7.61800 34.74300
C  40.24300 11.48300 34.00600
C  40.64100 12.04200 35.33400
H  40.55300 11.05000 36.56000
N  38.75400 7.55700 29.53000
C  38.00600 7.28200 30.63600
C  36.89200 6.50300 30.21600
C  36.97400 6.32300 28.78000
C  38.23000 7.03300 28.36500
C  35.87900 5.90900 31.21900
C  36.08100 5.65900 27.77900
O  36.28700 5.46300 26.57600
C  34.68000 5.19900 28.33400
N  40.98100 8.10600 27.39000
C  40.09800 7.49500 26.64000
C  40.46000 7.54300 25.15700
C  41.81600 8.30500 25.16100
C  42.00400 8.59400 26.64900
C  40.56900 6.08800 24.47000
C  41.66000 9.63200 24.29500
C  40.74000 10.78700 24.76200
N  42.48900 9.50900 29.49800
C  43.34000 9.72700 28.43100
C  44.44800 10.55700 28.94600
C  44.20300 10.86700 30.25000
C  42.96700 10.18300 30.54000
C  45.69600 10.84800 28.18600
C  44.66700 11.56800 31.44700
O  45.66700 12.23600 31.65000
C  43.53700 11.23400 32.51700
C  44.17900 10.46300 33.60100
O  44.82400 9.40400 33.50900
O  43.90700 11.19500 34.73900
C  44.28800 10.53200 35.99300
H  37.46200 7.40700 32.64000
H  38.23300 6.48700 26.39900
H  43.84700 9.59300 26.33600
H  41.53500 9.82500 34.41700
H  38.70200 9.27400 34.20100
H  40.95100 7.85100 35.37900
H  40.49700 6.86600 34.06300
H  39.17400 7.30300 35.22800
H  39.15600 11.54000 33.94800
H  40.73700 12.04000 33.21000
H  39.99600 12.85500 35.66500
H  41.60800 12.52200 35.48500
H  35.14900 6.60800 31.62700
H  36.53200 5.46300 31.96900
H  35.20300 5.18600 30.76300
H  34.82100 4.32000 28.96300
H  33.98600 5.15200 27.49500
H  34.17100 5.96300 28.92200
H  39.71300 8.14400 24.63900
H  42.64100 7.73100 24.73900
H  39.67900 5.98600 23.84800
H  40.55600 5.31000 25.23400
H  41.39700 5.86900 23.79700
H  41.39900 9.47200 23.24900
H  42.67000 10.03800 24.36000
H  41.45000 11.59500 24.93900
H  40.14900 10.52600 25.63900
H  40.01100 11.10100 24.01500
H  45.97000 11.83700 28.55200
H  45.54800 10.83600 27.10600
H  46.48500 10.23100 28.61400
H  43.32700 12.24400 32.86800
H  44.09900 11.31400 36.72900
H  45.36500 10.37600 35.93100
H  43.84100 9.55300 36.16700
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


