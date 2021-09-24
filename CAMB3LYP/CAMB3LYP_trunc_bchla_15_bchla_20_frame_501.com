%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 47.33900 34.91800 28.07800
C  45.68800 33.19800 30.68300
C  48.10100 37.44900 30.24400
C  48.24300 36.72400 25.47200
C  46.07200 32.47700 25.83700
N  47.02000 35.19800 30.25900
C  46.36000 34.38500 31.06000
C  46.67100 34.77400 32.49400
C  47.03300 36.30500 32.38100
C  47.41200 36.35800 30.85400
C  45.84200 37.21500 32.86800
C  47.88800 34.00400 33.06800
C  47.98300 33.94800 34.65600
H  46.95600 34.64300 35.55400
N  48.24300 36.70300 27.90700
C  48.56900 37.63000 28.93000
C  49.30300 38.76000 28.31400
C  49.26700 38.57700 26.88500
C  48.57000 37.29400 26.64900
C  49.93800 39.85800 29.09200
C  49.83100 39.42600 25.68000
O  49.78800 39.06500 24.49600
C  50.32700 40.81500 25.96000
N  47.07800 34.69300 25.91700
C  47.58500 35.54400 25.06500
C  47.43000 35.17200 23.61700
C  46.71000 33.78100 23.71200
C  46.68600 33.56800 25.25900
C  46.76800 36.29200 22.80900
C  47.43400 32.63900 23.01000
C  46.53900 31.60900 22.23200
N  46.16500 33.09500 28.16200
C  45.77400 32.19200 27.15400
C  45.18200 31.05300 27.76700
C  44.95200 31.44600 29.12400
C  45.59700 32.64100 29.31900
C  44.92100 29.73400 27.09800
C  44.48300 31.08700 30.48400
O  43.89600 30.09800 30.87100
C  44.73800 32.34900 31.43800
C  44.99200 31.90700 32.84200
O  45.89000 31.21500 33.22200
O  43.96300 32.32800 33.70500
C  43.67300 31.45600 34.86300
H  48.20800 38.37100 30.82000
H  48.63100 37.17100 24.55500
H  45.76300 31.70600 25.12800
H  45.78600 34.43600 33.03300
H  47.91000 36.55500 32.97900
H  44.99600 36.59900 33.17100
H  45.50600 37.70100 31.95200
H  46.07700 37.86100 33.71400
H  48.75200 34.55000 32.68700
H  47.86600 33.00400 32.63500
H  48.97800 34.33500 34.88100
H  48.01300 32.88400 34.89200
H  50.02200 39.76300 30.17400
H  49.40400 40.80200 28.98700
H  50.95300 39.89900 28.69600
H  51.37600 40.86200 26.25200
H  49.62900 41.09200 26.75100
H  50.11900 41.39300 25.06000
H  48.46800 35.04700 23.30900
H  45.69100 33.78300 23.32700
H  46.13600 36.95200 23.40400
H  46.11300 35.78300 22.10300
H  47.49100 36.91500 22.28200
H  48.02000 32.05200 23.71700
H  48.21400 33.13900 22.43600
H  45.47900 31.78000 22.41900
H  46.81500 30.64200 22.65200
H  46.65400 31.52300 21.15100
H  44.41500 29.99000 26.16700
H  44.28000 29.05700 27.66300
H  45.86600 29.26100 26.82900
H  43.79600 32.88200 31.56300
H  43.90900 32.01600 35.76800
H  44.22800 30.51800 34.82500
H  42.60300 31.25600 34.80700
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


