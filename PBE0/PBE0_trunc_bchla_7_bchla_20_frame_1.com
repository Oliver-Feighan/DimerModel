%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 26.07000 0.32300 29.17500
C  27.93000 0.00300 32.17200
C  23.23700 0.37300 31.22200
C  24.20000 0.73500 26.50800
C  28.89900 -0.55900 27.35600
N  25.68800 0.23300 31.42600
C  26.58700 0.25300 32.44900
C  25.99300 0.25100 33.86800
C  24.43500 0.24800 33.45400
C  24.48700 0.27000 31.95500
C  23.63400 -1.02300 33.94300
C  26.39900 1.50900 34.70300
C  27.30600 1.47800 35.92700
H  27.07500 2.41400 37.10600
N  24.00400 0.60600 28.87400
C  22.98600 0.52400 29.83500
C  21.72400 0.61300 29.13500
C  21.97300 0.88000 27.79100
C  23.44300 0.80400 27.68300
C  20.46800 0.64800 29.93700
C  21.02800 1.08200 26.60700
O  21.39200 1.28100 25.44400
C  19.58500 1.18400 26.91000
N  26.42900 0.02500 27.20000
C  25.57000 0.46200 26.23000
C  26.30100 0.60400 24.85200
C  27.70600 0.10100 25.15700
C  27.68000 -0.24400 26.67200
C  25.55900 -0.10800 23.67400
C  28.73300 1.18300 24.79000
C  29.96300 0.77300 24.06500
N  28.00500 -0.28200 29.61000
C  29.05600 -0.58000 28.78000
C  30.24500 -0.86800 29.52500
C  29.92400 -0.56000 30.93700
C  28.52000 -0.19800 30.89800
C  31.58600 -1.33700 28.93500
C  30.40700 -0.48700 32.33900
O  31.55800 -0.58600 32.79700
C  29.11500 -0.11100 33.17500
C  29.05800 -1.28400 34.10300
O  28.53100 -2.37500 33.92200
O  29.56000 -0.97800 35.26800
C  29.75400 -2.03800 36.28700
H  22.46700 0.52900 31.98000
H  23.69600 0.81600 25.54300
H  29.79300 -0.76400 26.76400
H  26.26300 -0.69600 34.33700
H  23.82700 1.11200 33.72400
H  24.21800 -1.88400 34.26700
H  22.99700 -1.39000 33.13900
H  22.99100 -0.72400 34.77100
H  25.46900 1.94600 35.06700
H  26.90100 2.23700 34.06700
H  28.23400 1.84700 35.49000
H  27.27900 0.42100 36.19200
H  19.92700 1.58400 29.80500
H  20.67800 0.53600 31.00100
H  19.76300 -0.13800 29.66800
H  19.43300 1.90200 27.71600
H  19.18900 0.27500 27.36300
H  18.97400 1.53500 26.07900
H  26.25300 1.68100 24.69200
H  27.96600 -0.74800 24.52500
H  24.56100 -0.46700 23.92500
H  26.08600 -0.97800 23.28400
H  25.56500 0.53500 22.79400
H  29.00400 1.58600 25.76500
H  28.34900 1.97700 24.15000
H  30.78000 1.14300 24.68400
H  29.93300 1.13400 23.03700
H  30.14200 -0.30200 24.03500
H  31.60000 -2.42600 28.91000
H  32.38900 -0.93100 29.55000
H  31.71700 -1.01400 27.90200
H  29.35100 0.80400 33.71800
H  28.83300 -2.53300 36.59300
H  30.27100 -1.70300 37.18600
H  30.41000 -2.86500 36.01700
Mg 7.14200 57.54600 41.51700
C  6.28900 54.20700 41.28700
C  10.28800 56.72100 40.38900
C  7.81200 60.89700 41.16900
C  3.94100 58.30400 42.53200
N  8.10000 55.73200 40.72500
C  7.55000 54.46800 40.80800
C  8.51900 53.39300 40.32300
C  9.88500 54.17400 40.37000
C  9.39700 55.61300 40.45200
C  10.94100 53.73800 41.48400
C  8.18100 52.71900 38.95500
C  8.31100 51.21100 38.91500
H  9.48100 50.75900 38.05500
N  8.80400 58.60500 40.85800
C  9.98200 58.10200 40.47200
C  10.93200 59.16600 40.19700
C  10.33500 60.40200 40.38100
C  8.90800 60.02900 40.77200
C  12.41400 58.88100 39.65800
C  11.01900 61.80700 40.27800
O  12.22600 61.86700 39.90900
C  10.28000 63.06100 40.55000
N  5.92900 59.42300 41.58100
C  6.46000 60.62400 41.57700
C  5.53200 61.67900 42.14200
C  4.20400 60.94400 42.52900
C  4.70000 59.43900 42.28400
C  6.08200 62.75400 43.13300
C  2.88100 61.24800 41.76800
C  1.64800 60.78800 42.51800
N  5.41300 56.46200 41.84200
C  4.16700 56.89800 42.34900
C  3.34900 55.73600 42.56600
C  4.13900 54.62000 42.19700
C  5.33900 55.15600 41.73800
C  2.01600 55.69000 43.29800
C  4.11700 53.24300 42.14700
O  3.28800 52.39200 42.32200
C  5.50100 52.88200 41.39400
C  5.07100 52.59100 39.95700
O  4.83200 53.42500 39.10500
O  4.99200 51.21300 39.83600
C  4.48800 50.71300 38.56100
H  11.34700 56.50800 40.23000
H  8.09700 61.94400 41.04900
H  2.94700 58.37800 42.97800
H  8.57300 52.63000 41.09900
H  10.29300 54.14500 39.35900
H  11.53300 54.64700 41.59000
H  11.67500 53.00700 41.14400
H  10.53800 53.39000 42.43500
H  8.89800 53.17500 38.27200
H  7.26400 53.03800 38.46000
H  7.37100 50.81600 38.52900
H  8.48000 50.78000 39.90200
H  12.65000 57.83100 39.48000
H  13.15600 59.27900 40.35100
H  12.59000 59.28000 38.66000
H  9.90800 63.04900 41.57400
H  9.55300 63.27700 39.76700
H  11.06400 63.81300 40.63000
H  5.09100 62.35700 41.41100
H  4.06100 60.88100 43.60800
H  7.16900 62.67900 43.08000
H  5.65900 62.68700 44.13600
H  5.85900 63.78200 42.84700
H  2.99200 60.79700 40.78200
H  2.91000 62.32500 41.60300
H  0.83700 61.47100 42.26400
H  1.73700 60.78400 43.60500
H  1.32700 59.82700 42.11600
H  1.67100 54.68100 43.52400
H  1.28900 56.22800 42.69100
H  2.06800 56.22200 44.24800
H  6.04900 52.10100 41.92200
H  3.44600 51.01000 38.44100
H  4.61600 49.63500 38.46900
H  5.02300 51.19800 37.74500


