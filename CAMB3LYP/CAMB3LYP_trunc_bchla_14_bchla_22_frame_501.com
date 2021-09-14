%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 46.48800 44.18100 43.47500
C  43.15300 43.46100 43.30400
C  47.10800 40.75600 42.73500
C  49.77100 44.92500 43.24200
C  45.73500 47.57600 43.26500
N  45.18800 42.33600 42.67700
C  43.85400 42.27500 42.94800
C  43.32200 40.91600 42.58600
C  44.62300 40.02800 42.36600
C  45.76600 41.06300 42.61200
C  44.58600 38.65100 42.92100
C  42.20900 40.77400 41.44300
C  42.22600 41.75800 40.31100
H  42.13200 41.19200 38.90100
N  48.19100 42.96300 43.28600
C  48.25100 41.62900 43.02600
C  49.60100 41.16100 42.91600
C  50.41700 42.37800 43.18500
C  49.45800 43.50200 43.25700
C  50.07700 39.81200 42.65600
C  51.92500 42.49800 43.35200
O  52.62900 41.58000 43.12100
C  52.68900 43.74400 43.70500
N  47.56100 45.95400 43.15400
C  48.89600 45.98900 43.28400
C  49.46200 47.47200 42.97200
C  48.15100 48.31100 43.36800
C  47.10700 47.23300 43.19900
C  50.71200 47.88700 43.75300
C  48.00000 49.60500 42.54800
C  47.68500 50.93800 43.23700
N  44.79500 45.31200 43.20500
C  44.59800 46.68700 43.17600
C  43.22600 47.01600 43.20700
C  42.57500 45.79800 43.31300
C  43.57500 44.76200 43.30900
C  42.60800 48.38000 43.21500
C  41.37100 45.09400 43.47000
O  40.23700 45.56500 43.52900
C  41.64700 43.57300 43.29900
C  40.96100 42.96000 44.51100
O  40.13800 42.02200 44.45100
O  41.52200 43.52500 45.63200
C  41.27100 42.99100 46.99300
H  47.22800 39.72900 42.38400
H  50.80000 45.16300 42.96500
H  45.50900 48.64400 43.27100
H  42.74600 40.53700 43.43100
H  44.70600 39.88500 41.28900
H  45.47200 38.39300 43.50100
H  44.34100 37.90700 42.16300
H  43.78400 38.58200 43.65600
H  41.18300 40.86000 41.80100
H  42.35800 39.83500 40.90900
H  42.98200 42.53800 40.39500
H  41.27000 42.28200 40.33700
H  49.35900 39.10300 43.06700
H  51.06100 39.59900 43.07200
H  50.12200 39.70100 41.57200
H  51.89900 44.16700 44.32600
H  52.89000 44.37800 42.84200
H  53.62800 43.61700 44.24400
H  49.58000 47.51900 41.88900
H  48.02200 48.62000 44.40500
H  50.38700 48.59700 44.51400
H  51.51000 48.17200 43.06800
H  51.16800 47.06300 44.30100
H  47.21200 49.41800 41.81800
H  48.88600 49.70600 41.92100
H  48.29300 51.80900 42.99200
H  47.59300 50.84600 44.32000
H  46.63600 51.07800 42.97800
H  41.68900 48.46100 43.79600
H  42.39500 48.66700 42.18500
H  43.24900 49.14100 43.66200
H  41.15500 43.34600 42.35400
H  42.16000 43.01500 47.62300
H  40.93100 41.96000 46.89400
H  40.51300 43.43800 47.63500
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


