%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 3.63500 0.41800 44.38500
C  6.38200 2.41600 44.15700
C  1.85600 2.91600 42.76500
C  1.03400 -1.70200 44.27100
C  5.65900 -2.27700 45.39800
N  4.15300 2.34100 43.32400
C  5.35700 2.93700 43.43600
C  5.44100 4.20100 42.66100
C  3.82600 4.53100 42.53700
C  3.19500 3.17000 42.89900
C  3.38200 5.72400 43.43400
C  6.08500 3.93300 41.31000
C  5.45900 2.81500 40.35400
H  5.36900 3.04800 38.85100
N  1.72200 0.57100 43.63100
C  1.11700 1.74400 43.11300
C  -0.27500 1.47100 42.97900
C  -0.55500 0.16000 43.34400
C  0.72000 -0.40500 43.81700
C  -1.17200 2.58100 42.47000
C  -1.87000 -0.57100 43.21500
O  -2.81800 0.09700 42.73200
C  -2.07300 -1.97000 43.39000
N  3.45400 -1.75300 44.57800
C  2.23500 -2.36400 44.54600
C  2.35900 -3.87900 44.83600
C  3.87500 -4.03200 45.19700
C  4.39100 -2.62000 44.97700
C  1.41200 -4.36800 45.91100
C  4.59000 -5.14900 44.26500
C  5.56000 -6.06500 45.02000
N  5.66900 0.04500 44.75000
C  6.32000 -0.99300 45.33600
C  7.63200 -0.52000 45.68500
C  7.75100 0.85800 45.29200
C  6.52700 1.13100 44.70400
C  8.77000 -1.38000 46.11700
C  8.57600 2.01500 45.15500
O  9.73800 2.10000 45.57200
C  7.71800 3.05600 44.46000
C  7.61900 4.24800 45.33700
O  8.63700 4.80700 45.69500
O  6.34700 4.66100 45.59200
C  6.23700 5.85100 46.43800
H  1.24500 3.78600 42.51600
H  0.13600 -2.29700 44.44900
H  6.31200 -3.00800 45.88000
H  6.00100 5.03200 43.09000
H  3.58600 4.93100 41.55100
H  2.93500 6.50500 42.82000
H  4.28700 6.13300 43.88400
H  2.70200 5.45700 44.24400
H  7.14700 3.76700 41.49200
H  6.12500 4.93200 40.87600
H  4.43400 2.52000 40.57900
H  6.17000 1.99900 40.48000
H  -1.10900 2.52700 41.38300
H  -0.76100 3.53900 42.78600
H  -2.21500 2.54500 42.78600
H  -2.88500 -2.34300 42.76500
H  -2.15000 -2.27300 44.43500
H  -1.20000 -2.43300 42.93200
H  2.14200 -4.33200 43.86800
H  3.97500 -4.21900 46.26600
H  0.70500 -3.61400 46.25600
H  1.77700 -4.91000 46.78300
H  0.73700 -5.01500 45.35000
H  5.15900 -4.73200 43.43400
H  3.86400 -5.76700 43.73700
H  5.09600 -7.04700 45.10800
H  5.65500 -5.69200 46.04000
H  6.52700 -6.16100 44.52600
H  8.71100 -1.83500 47.10500
H  9.65100 -0.74100 46.04600
H  8.91100 -2.23600 45.45700
H  8.25000 3.32500 43.54700
H  6.78900 6.65600 45.95200
H  6.55900 5.68400 47.46600
H  5.19100 6.09800 46.61900
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


