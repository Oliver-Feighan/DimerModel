%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 44.01800 2.91800 47.66100
C  42.06000 5.71500 47.01900
C  41.00600 1.01600 47.31200
C  45.79900 0.15700 47.76200
C  46.87900 4.91700 47.24100
N  41.83200 3.33000 47.09800
C  41.29400 4.59200 46.98100
C  39.79500 4.50000 46.68300
C  39.48700 3.02600 47.06600
C  40.84700 2.38200 47.12900
C  38.69200 2.74400 48.40100
C  39.44800 4.83900 45.22400
C  38.00300 4.93400 44.91400
H  37.59500 5.76800 43.64000
N  43.49300 0.83400 47.54100
C  42.21400 0.26700 47.42800
C  42.35900 -1.16200 47.43400
C  43.77600 -1.43800 47.42700
C  44.43200 -0.09700 47.60800
C  41.15600 -2.13400 47.28600
C  44.34200 -2.85900 47.34100
O  43.62700 -3.84400 47.42200
C  45.80100 -3.04100 47.16300
N  46.06200 2.62900 47.74100
C  46.53800 1.38200 47.82400
C  48.11500 1.40300 47.69600
C  48.51000 2.91300 47.47000
C  47.02000 3.56100 47.42600
C  48.93400 0.68700 48.80200
C  49.24600 3.07700 46.12300
C  50.25900 4.22700 46.06800
N  44.47400 4.92000 47.22100
C  45.65200 5.61300 47.20700
C  45.34700 7.04400 47.09000
C  43.98000 7.11700 47.10900
C  43.47900 5.83500 47.12000
C  46.36600 8.18800 47.10600
C  42.87400 8.02300 47.02800
O  42.81600 9.22700 46.97900
C  41.57700 7.14800 46.98100
C  40.79100 7.44500 48.21000
O  41.19100 7.46200 49.34100
O  39.43700 7.57100 47.85400
C  38.52600 7.78100 48.94300
H  40.02600 0.53600 47.28500
H  46.51600 -0.66500 47.82200
H  47.74600 5.57800 47.17700
H  39.30800 5.14400 47.41400
H  38.96600 2.43700 46.31200
H  37.88000 2.08800 48.08800
H  38.53600 3.57700 49.08700
H  39.33000 2.09900 49.00500
H  39.75400 4.10000 44.48300
H  39.93900 5.76800 44.93400
H  37.52300 5.35800 45.79600
H  37.59600 3.92700 44.82400
H  41.34100 -2.41600 46.24900
H  40.16900 -1.67200 47.25600
H  41.01800 -3.02000 47.90600
H  46.05900 -2.38200 46.33400
H  45.98200 -4.06600 46.83900
H  46.21200 -2.73700 48.12600
H  48.33800 0.83400 46.79300
H  49.07900 3.37000 48.27900
H  49.32800 -0.23900 48.38200
H  48.30000 0.56800 49.68000
H  49.82500 1.27000 49.03500
H  48.45300 3.22300 45.39000
H  49.73200 2.16800 45.76800
H  51.18600 3.72000 45.79800
H  50.32400 4.81900 46.98100
H  50.02100 4.90100 45.24500
H  46.22000 8.97700 46.36900
H  47.30400 7.65100 46.96300
H  46.44100 8.50600 48.14600
H  41.03900 7.44900 46.08200
H  38.18200 8.81200 48.86000
H  38.86300 7.60600 49.96500
H  37.69400 7.09100 48.80500
Mg 7.63600 56.80200 41.94500
C  6.62500 53.44600 41.38500
C  10.68100 55.97500 40.45800
C  8.03200 60.02500 41.40500
C  4.22800 57.38500 43.14900
N  8.45000 54.91300 40.78100
C  7.95700 53.66100 40.95000
C  8.99000 52.64800 40.49000
C  10.32900 53.33300 40.76400
C  9.81400 54.85400 40.60800
C  10.97600 53.14900 42.13800
C  8.68800 52.22300 38.98700
C  8.78600 50.76200 38.78500
H  9.47200 50.26400 37.49800
N  9.15000 57.91900 40.99500
C  10.37200 57.38400 40.51700
C  11.15700 58.50500 40.06500
C  10.49600 59.72000 40.38700
C  9.15900 59.29200 40.93300
C  12.51600 58.28800 39.40700
C  10.98000 61.08400 40.08900
O  12.05800 61.21700 39.52300
C  10.25300 62.36800 40.32000
N  6.19300 58.49600 42.17300
C  6.76600 59.71000 42.05200
C  5.81300 60.81300 42.54700
C  4.61500 59.97500 43.27000
C  5.02400 58.52600 42.89900
C  6.51300 61.95400 43.41400
C  3.16600 60.37300 42.84800
C  2.45300 59.49100 41.76700
N  5.85400 55.73100 42.22300
C  4.59600 56.10000 42.77900
C  3.77800 54.83100 43.01200
C  4.55200 53.84100 42.42800
C  5.79000 54.43400 41.99000
C  2.39800 54.58300 43.65100
C  4.62700 52.43800 42.17000
O  3.86200 51.52600 42.48800
C  5.87100 52.24300 41.37900
C  5.33600 51.93400 39.96200
O  4.91100 52.73800 39.13500
O  5.35200 50.52000 39.80700
C  4.79600 49.97200 38.63700
H  11.68700 55.76400 40.08900
H  8.33200 61.07300 41.34400
H  3.25200 57.43300 43.63500
H  9.01800 51.75000 41.10700
H  11.18700 53.14700 40.11900
H  10.60100 52.28700 42.68900
H  11.03000 54.09000 42.68600
H  12.00300 52.89400 41.87400
H  9.43500 52.68300 38.34000
H  7.77500 52.65300 38.57500
H  7.73900 50.46500 38.83400
H  9.33600 50.33800 39.62500
H  13.35000 58.53400 40.06400
H  12.53800 58.83400 38.46500
H  12.64500 57.23800 39.14400
H  10.18000 62.52800 41.39600
H  9.25600 62.43600 39.88300
H  10.86800 63.22900 40.06100
H  5.41300 61.22100 41.61900
H  4.75100 60.04000 44.34900
H  7.57600 62.02500 43.18400
H  6.34900 61.97700 44.49200
H  5.91600 62.83300 43.16600
H  3.09900 61.41000 42.52000
H  2.61700 60.24800 43.78100
H  1.86900 60.30400 41.33500
H  1.78000 58.70400 42.10800
H  3.21700 59.14200 41.07200
H  2.47400 53.89100 44.48900
H  1.98400 53.99700 42.83100
H  1.91500 55.55200 43.77500
H  6.40600 51.37000 41.75300
H  4.82800 48.88300 38.58500
H  5.36500 50.36200 37.79400
H  3.76000 50.31200 38.63400


