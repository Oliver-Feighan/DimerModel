%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 40.80400 41.48400 27.15800
C  39.93800 43.96600 29.51800
C  41.67100 39.43500 29.59200
C  42.41200 39.59100 24.81200
C  40.80000 44.20800 24.65700
N  40.84000 41.74700 29.31900
C  40.23400 42.70800 30.04600
C  40.11400 42.23500 31.50800
C  40.62500 40.83000 31.52200
C  41.12800 40.63200 30.06000
C  41.62900 40.55000 32.71400
C  38.74000 42.40700 32.23800
C  38.20000 41.35400 33.29800
H  37.58800 41.87900 34.60500
N  41.67300 39.62800 27.19200
C  41.94100 38.92200 28.34300
C  42.50600 37.66000 27.95800
C  42.79100 37.65600 26.53700
C  42.34000 39.00600 26.13300
C  42.84100 36.69600 29.02300
C  43.34600 36.53900 25.61900
O  43.54900 36.75600 24.47000
C  43.39500 35.12400 26.09600
N  41.54800 41.90500 25.04000
C  42.04300 40.86300 24.33500
C  42.34300 41.24800 22.85100
C  41.79900 42.70700 22.76800
C  41.36100 43.00200 24.22400
C  43.78600 41.15300 22.43200
C  40.63300 43.00600 21.76200
C  39.24100 42.61000 22.12200
N  40.45200 43.63600 26.97900
C  40.44000 44.51600 25.93500
C  39.94600 45.79200 26.39200
C  39.72700 45.66500 27.78900
C  40.09700 44.34100 28.10100
C  39.86300 47.06500 25.62600
C  39.33400 46.33800 28.99400
O  38.92600 47.46900 29.23200
C  39.40400 45.25700 30.16400
C  40.16300 45.75100 31.38400
O  41.33900 46.00900 31.38400
O  39.30000 45.80400 32.46700
C  39.93400 46.30500 33.68700
H  41.79100 38.62500 30.31500
H  42.89000 38.90600 24.10900
H  40.69600 44.98800 23.90000
H  40.78700 42.96300 31.96100
H  39.78000 40.17800 31.74500
H  41.64600 41.41100 33.38200
H  42.57700 40.23500 32.27800
H  41.25100 39.76900 33.37400
H  37.98000 42.45900 31.45800
H  38.66600 43.37100 32.74200
H  39.06500 40.72800 33.51800
H  37.46000 40.71500 32.81600
H  43.75800 36.20700 28.69500
H  42.05100 35.94600 29.03300
H  42.98100 37.00500 30.05900
H  44.39900 35.05700 26.51300
H  43.13900 34.46100 25.26900
H  42.74200 34.89700 26.93900
H  41.71500 40.69400 22.15300
H  42.49200 43.53200 22.60600
H  44.38600 40.63600 23.18200
H  44.24800 42.13400 22.31800
H  43.78100 40.59500 21.49600
H  40.84900 42.53700 20.80200
H  40.70100 44.05800 21.48400
H  39.26600 41.80800 22.86000
H  38.56300 42.19300 21.37700
H  38.76600 43.53300 22.45600
H  40.81300 47.42600 25.23300
H  39.70300 47.84000 26.37600
H  39.09100 47.02600 24.85700
H  38.36800 45.04000 30.42500
H  39.84500 45.31900 34.14200
H  39.21300 47.03100 34.06200
H  40.94500 46.67900 33.52400
Mg 9.20000 47.36600 24.88000
C  7.02600 47.76900 27.74700
C  11.57000 48.44700 26.79100
C  10.94600 47.19300 22.09000
C  6.29100 46.96700 22.73800
N  9.27600 47.92100 27.10800
C  8.32000 48.03400 28.09600
C  8.93500 48.52700 29.35600
C  10.25500 49.18700 28.88100
C  10.38000 48.52600 27.48800
C  10.19000 50.72500 28.75600
C  9.27100 47.32000 30.28400
C  8.91000 47.48600 31.81100
H  8.76300 46.21000 32.73200
N  11.12300 47.75000 24.52700
C  11.99300 48.08000 25.44400
C  13.30700 48.11300 24.83000
C  13.09300 47.62600 23.51900
C  11.67000 47.53200 23.33600
C  14.56100 48.44000 25.66100
C  14.15200 47.33300 22.49300
O  13.84400 46.95300 21.40300
C  15.59100 47.42400 22.88400
N  8.59400 47.38800 22.62900
C  9.59400 47.15300 21.78200
C  9.10700 46.72400 20.44900
C  7.56400 46.85800 20.64800
C  7.41500 47.09900 22.12700
C  9.70100 47.44000 19.18300
C  6.80800 45.57800 20.11700
C  6.00900 45.65200 18.78600
N  7.17400 47.35700 25.08600
C  6.10400 47.11400 24.19600
C  4.83100 47.11900 24.94600
C  5.17200 47.34400 26.27300
C  6.54300 47.52700 26.33600
C  3.42100 47.05000 24.31100
C  4.55900 47.48400 27.63300
O  3.43200 47.38700 28.05800
C  5.73500 47.94200 28.63300
C  5.68600 47.23700 29.86100
O  6.32700 46.19300 30.04600
O  4.84600 47.74900 30.84500
C  4.90400 47.16700 32.22100
H  12.46500 48.75400 27.33500
H  11.53100 47.07200 21.17600
H  5.40800 46.83900 22.11000
H  8.32300 49.32300 29.78100
H  11.13800 49.00400 29.49300
H  10.92100 51.19300 29.41600
H  9.18200 51.02300 29.04300
H  10.33700 51.10600 27.74500
H  10.35100 47.19900 30.19800
H  8.78700 46.45500 29.83100
H  8.00100 48.08400 31.87000
H  9.74500 48.08400 32.17800
H  15.18900 47.54900 25.70000
H  14.38200 48.83100 26.66200
H  15.20000 49.13500 25.11700
H  15.85600 46.71400 23.66700
H  15.91500 48.35300 23.35300
H  16.20400 47.15000 22.02500
H  9.34100 45.66500 20.34000
H  7.22100 47.78300 20.18400
H  10.49600 48.09000 19.54700
H  9.00000 47.94700 18.51900
H  10.12400 46.67900 18.52700
H  6.07000 45.27200 20.85900
H  7.53000 44.76900 20.22800
H  6.23700 44.78300 18.16900
H  6.25500 46.57800 18.26600
H  4.93600 45.62800 18.98000
H  3.49200 46.43400 23.41500
H  3.14200 48.04200 23.95500
H  2.68000 46.64400 25.00000
H  5.61100 49.00700 28.82700
H  5.88500 47.32300 32.67200
H  4.47300 46.16700 32.25100
H  4.32200 47.83100 32.86000


