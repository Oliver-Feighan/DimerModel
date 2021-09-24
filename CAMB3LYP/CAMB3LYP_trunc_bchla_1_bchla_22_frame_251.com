%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg -2.04400 17.49100 26.73800
C  -2.15500 15.53600 29.75200
C  -2.48400 20.32400 28.82400
C  -2.10000 19.56800 23.93100
C  -2.02900 14.73500 24.85000
N  -2.50300 17.81300 29.09300
C  -2.37100 16.88400 30.06400
C  -2.66300 17.53000 31.46600
C  -2.96600 19.01100 31.04000
C  -2.55400 19.10600 29.55300
C  -4.45800 19.41500 31.30700
C  -1.49000 17.42200 32.44000
C  -1.82900 17.03100 33.84100
H  -0.68700 16.80800 34.84300
N  -2.05000 19.60600 26.40500
C  -2.24700 20.56500 27.39300
C  -2.18900 21.87200 26.76900
C  -1.98200 21.71100 25.40900
C  -2.11500 20.23400 25.17200
C  -2.46700 23.16900 27.47300
C  -1.85900 22.77200 24.24000
O  -1.68000 22.55800 23.06600
C  -1.78900 24.24000 24.60100
N  -2.21300 17.23900 24.73600
C  -2.21800 18.17300 23.73300
C  -2.50700 17.54400 22.35000
C  -2.23600 16.00100 22.60900
C  -2.12000 15.99700 24.18100
C  -3.92800 17.89500 21.80400
C  -0.96100 15.44100 21.90800
C  0.37100 15.96600 22.43900
N  -2.11100 15.55600 27.12300
C  -2.02800 14.45800 26.22200
C  -2.06000 13.19200 26.95100
C  -2.05400 13.57600 28.30400
C  -2.07200 14.99400 28.37800
C  -2.02800 11.85900 26.36100
C  -2.05700 13.00500 29.64700
O  -1.98500 11.83800 29.99400
C  -2.21400 14.25700 30.66700
C  -1.01500 13.99300 31.58100
O  0.19500 14.13800 31.33000
O  -1.45400 13.50000 32.72600
C  -0.50700 13.41700 33.91100
H  -2.54500 21.19300 29.48200
H  -2.11900 20.13500 22.99700
H  -1.94800 13.88000 24.17600
H  -3.47200 16.92500 31.87400
H  -2.32200 19.72000 31.56200
H  -4.61000 20.30200 31.92200
H  -5.09400 18.58600 31.61800
H  -5.01600 19.78600 30.44700
H  -1.29200 18.47500 32.64100
H  -0.61800 16.92800 32.01000
H  -2.49700 16.18800 34.02500
H  -2.37000 17.90200 34.21100
H  -2.54100 23.02800 28.55100
H  -3.37300 23.68600 27.15300
H  -1.69300 23.93100 27.38000
H  -0.93800 24.47000 25.24200
H  -2.76500 24.39500 25.06100
H  -1.83600 24.73500 23.63100
H  -1.82000 18.06600 21.68300
H  -3.08800 15.32400 22.54300
H  -3.89400 18.55600 20.93800
H  -4.67200 18.26800 22.50700
H  -4.29000 16.90500 21.52600
H  -1.14700 15.45600 20.83400
H  -0.92500 14.37100 22.11200
H  0.30800 16.76300 23.17900
H  0.88600 16.37900 21.57200
H  0.96700 15.25500 23.01000
H  -2.06800 11.03700 27.07600
H  -1.09600 11.82400 25.79700
H  -2.85300 11.79600 25.65000
H  -3.07800 14.15100 31.32300
H  -1.14400 13.30000 34.78800
H  0.04800 14.34000 34.07600
H  0.23900 12.62700 33.81700
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


