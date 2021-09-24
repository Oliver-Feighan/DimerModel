%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
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
Mg -2.51200 34.32000 26.71000
C  -3.70600 32.37000 29.49400
C  -1.09600 36.45200 28.85600
C  -2.23500 36.47200 24.15600
C  -4.45400 32.24200 24.65700
N  -2.47200 34.41100 28.88300
C  -2.84200 33.43600 29.82500
C  -2.27900 33.77800 31.25700
C  -1.51600 35.11400 31.00700
C  -1.64900 35.32300 29.47500
C  -1.96800 36.32500 31.78300
C  -1.39500 32.57400 31.88800
C  -1.81300 32.12300 33.29300
H  -0.66900 32.23700 34.36800
N  -1.63700 36.15800 26.50500
C  -1.11000 36.89500 27.48800
C  -0.61000 38.17800 26.99100
C  -1.04200 38.20400 25.62000
C  -1.69000 36.91400 25.36100
C  0.10900 39.14600 27.82900
C  -0.90500 39.45100 24.69100
O  -1.22300 39.41800 23.51400
C  -0.35500 40.75800 25.16600
N  -3.22300 34.33100 24.76300
C  -2.85500 35.30500 23.84100
C  -3.39300 34.95800 22.41400
C  -3.78700 33.44000 22.55100
C  -3.77600 33.28900 24.07200
C  -4.45800 35.91300 21.78300
C  -2.78800 32.46000 21.87700
C  -3.33300 31.25600 21.15100
N  -3.86100 32.67700 26.98500
C  -4.55900 31.95500 26.05900
C  -5.21600 30.76800 26.76200
C  -4.86600 30.86600 28.08700
C  -4.10200 32.07900 28.18800
C  -5.97600 29.69500 26.01600
C  -4.98700 30.30300 29.42000
O  -5.62100 29.26000 29.68400
C  -4.28200 31.35300 30.46900
C  -5.23900 32.03400 31.42500
O  -6.23700 32.64900 31.00000
O  -4.96700 31.73200 32.76000
C  -6.09800 32.02100 33.63400
H  -0.51900 37.13000 29.48800
H  -2.04400 37.20300 23.36800
H  -4.93400 31.57200 23.94000
H  -3.19000 33.89300 31.84400
H  -0.52800 34.94900 31.43900
H  -2.83400 36.13100 32.41600
H  -2.36000 36.92300 30.96100
H  -1.18900 36.74500 32.42000
H  -0.38600 32.98500 31.87000
H  -1.36700 31.74700 31.17900
H  -2.21200 31.11400 33.18900
H  -2.62300 32.76200 33.64700
H  -0.25500 40.16900 27.91900
H  1.01900 39.36100 27.27000
H  0.49900 38.75500 28.76900
H  0.66500 40.53600 25.48200
H  -0.91700 41.02300 26.06200
H  -0.46600 41.44500 24.32800
H  -2.50600 35.02500 21.78400
H  -4.81400 33.23900 22.24700
H  -4.50900 36.75000 22.47900
H  -5.42000 35.44200 21.58300
H  -4.08600 36.40600 20.88400
H  -2.24500 32.16900 22.77600
H  -2.06500 32.93000 21.21000
H  -2.75000 31.18100 20.23300
H  -4.37700 31.36900 20.85700
H  -3.03100 30.39100 21.74200
H  -5.97600 29.97600 24.96300
H  -7.01700 29.73400 26.33600
H  -5.54700 28.71900 26.24200
H  -3.48200 30.82800 30.99100
H  -6.59400 31.14900 34.05800
H  -6.86800 32.64200 33.17600
H  -5.69700 32.60500 34.46200


