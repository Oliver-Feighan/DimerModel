%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 46.36700 34.87300 28.25800
C  44.78600 32.87200 30.57800
C  46.38600 37.42000 30.77600
C  47.55600 36.98300 26.07000
C  45.97900 32.44300 25.82700
N  45.44200 35.14400 30.41100
C  45.07300 34.11700 31.15800
C  45.12500 34.48000 32.66900
C  45.08600 36.04500 32.56800
C  45.67500 36.28100 31.16800
C  43.63100 36.61800 32.70900
C  46.46000 33.91200 33.20000
C  46.69200 34.00900 34.71400
H  45.77100 34.85100 35.57400
N  47.01100 36.84800 28.33900
C  46.96000 37.70200 29.49400
C  47.69800 38.90500 29.13400
C  48.14000 38.73300 27.78500
C  47.54200 37.48800 27.30500
C  47.97900 40.11100 30.10600
C  49.09900 39.67100 26.98000
O  49.27000 39.43700 25.81500
C  49.74700 40.84400 27.64100
N  46.80600 34.69600 26.23000
C  47.30800 35.72100 25.50500
C  47.46400 35.39100 24.02100
C  47.08600 33.89200 23.98700
C  46.58100 33.65400 25.46800
C  46.60500 36.26500 23.03500
C  48.11600 32.92900 23.45400
C  49.59100 32.99000 23.91900
N  45.55600 33.03000 28.12500
C  45.46900 32.14000 27.09300
C  44.91700 30.91500 27.63200
C  44.64100 31.10300 29.03000
C  45.05100 32.44500 29.25900
C  44.67700 29.67800 26.82400
C  44.09100 30.54400 30.21900
O  43.56800 29.43800 30.47700
C  43.93800 31.83400 31.19500
C  44.25300 31.34600 32.61000
O  45.22200 30.77300 32.93400
O  43.18900 31.65700 33.40800
C  43.38400 31.62900 34.84300
H  46.32800 38.21100 31.52600
H  47.90200 37.66100 25.28700
H  45.81100 31.64300 25.10400
H  44.23500 34.14900 33.20400
H  45.76700 36.36500 33.35700
H  43.51000 37.05600 33.69900
H  42.93300 35.79700 32.54400
H  43.41600 37.43400 32.01900
H  47.26100 34.50400 32.75800
H  46.66000 32.91100 32.81900
H  47.74400 34.11200 34.98100
H  46.57400 32.98400 35.06600
H  47.70900 41.07100 29.66500
H  49.03500 40.03600 30.36500
H  47.45400 39.94800 31.04800
H  50.66300 40.94800 27.05900
H  50.08200 40.62400 28.65400
H  49.16800 41.73600 27.88000
H  48.52200 35.54500 23.81000
H  46.16900 33.69700 23.43100
H  46.31700 35.59300 22.22700
H  47.32700 37.01000 22.70000
H  45.80600 36.77400 23.57400
H  48.21500 33.30700 22.43700
H  47.86400 31.86800 23.45400
H  50.35800 33.20100 23.17500
H  49.94100 32.01500 24.25800
H  49.69200 33.75700 24.68700
H  44.20400 28.97200 27.50600
H  45.62900 29.23800 26.52400
H  44.15300 29.83100 25.88100
H  42.91200 32.17200 31.04800
H  42.52900 30.96500 34.97100
H  43.38800 32.61600 35.30600
H  44.28100 31.10000 35.16600
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


