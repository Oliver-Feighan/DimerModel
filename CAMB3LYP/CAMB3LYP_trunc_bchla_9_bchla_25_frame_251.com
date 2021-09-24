%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 35.66000 1.74600 29.23200
C  33.47200 2.34100 31.89500
C  38.12100 1.24100 31.62500
C  37.79000 1.61500 26.77900
C  33.01700 2.53400 27.01000
N  35.87600 2.18500 31.54100
C  34.82500 2.24700 32.36500
C  35.34000 2.18200 33.76400
C  36.88000 2.27300 33.63100
C  37.00600 1.84600 32.18200
C  37.59500 3.59000 33.93800
C  34.89600 0.87200 34.48900
C  34.26200 0.97200 35.91400
H  34.68800 2.14100 36.85600
N  37.69800 1.31100 29.19900
C  38.45500 0.95400 30.30700
C  39.74800 0.60100 29.81800
C  39.73600 0.82500 28.38300
C  38.39000 1.25600 28.03600
C  40.94400 0.20000 30.70900
C  40.89400 0.72300 27.34900
O  40.78100 0.94900 26.18600
C  42.22500 0.23800 27.88900
N  35.39600 2.02600 27.20700
C  36.47200 1.98000 26.39000
C  36.13800 2.22700 24.89300
C  34.55100 2.25100 24.96900
C  34.26500 2.24500 26.46600
C  36.85300 3.59900 24.36100
C  33.86400 1.14100 24.11300
C  33.41700 -0.12600 24.82000
N  33.60700 2.29400 29.34300
C  32.66300 2.48800 28.37100
C  31.34700 2.76700 29.04600
C  31.64200 2.71200 30.42700
C  33.04500 2.46100 30.58700
C  30.05700 3.20300 28.40600
C  31.07900 2.83900 31.76900
O  29.91300 3.09600 32.08200
C  32.20100 2.44500 32.72500
C  32.26400 3.54000 33.75100
O  32.77700 4.63100 33.58200
O  31.49900 3.13900 34.75300
C  31.03100 4.13600 35.70300
H  38.92600 1.22300 32.36200
H  38.44500 1.39800 25.93300
H  32.19900 2.76200 26.32300
H  34.93300 3.02400 34.32300
H  37.28500 1.57100 34.36000
H  36.96900 4.41100 34.28800
H  38.13600 3.93900 33.05800
H  38.32400 3.36100 34.71500
H  35.66700 0.10200 34.50700
H  34.14500 0.53000 33.77700
H  34.54900 0.03000 36.38300
H  33.17800 1.00900 35.80900
H  41.84000 0.79100 30.51900
H  41.29800 -0.80900 30.49700
H  40.73900 0.22800 31.78000
H  42.13800 -0.68600 28.46100
H  42.66900 1.05000 28.46500
H  42.94900 0.14100 27.08000
H  36.48000 1.41700 24.24900
H  34.11700 3.15500 24.54200
H  36.02900 4.30300 24.24400
H  37.25500 3.51600 23.35100
H  37.53000 4.05700 25.08300
H  34.65800 0.79900 23.44900
H  33.05000 1.52900 23.50100
H  34.25200 -0.36700 25.47900
H  33.30700 -0.80900 23.97800
H  32.56000 -0.11800 25.49300
H  30.04900 3.41600 27.33700
H  29.71100 4.11800 28.88700
H  29.33300 2.43700 28.68200
H  31.98900 1.57300 33.34400
H  31.22900 5.19300 35.52800
H  31.60100 3.89200 36.60000
H  30.07700 3.79700 36.10700
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


