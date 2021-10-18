%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 1.29600 7.71300 26.88100
C  1.46900 9.82700 29.67000
C  2.42600 5.06900 28.94700
C  1.47300 5.72500 24.13400
C  0.92500 10.50000 24.87700
N  1.69900 7.49100 29.04900
C  1.51500 8.49600 30.00300
C  1.43800 7.84300 31.43000
C  1.83300 6.38800 31.15200
C  2.05200 6.30500 29.60100
C  3.12000 6.01700 31.93700
C  0.04000 8.06800 31.99300
C  -0.12000 8.13500 33.50300
H  1.05800 7.66300 34.38500
N  1.84800 5.68100 26.53500
C  2.31100 4.77600 27.50900
C  2.63500 3.49600 26.90400
C  2.31300 3.59100 25.49300
C  1.87500 5.04400 25.33500
C  3.03900 2.24400 27.69300
C  2.37500 2.55400 24.38800
O  2.03400 2.84400 23.24900
C  2.64700 1.14500 24.75200
N  1.13700 8.09300 24.78900
C  1.15400 7.08500 23.87300
C  0.77800 7.64200 22.47500
C  0.53200 9.18300 22.70900
C  0.76500 9.29400 24.23800
C  1.89300 7.33600 21.45900
C  -0.89600 9.67800 22.35700
C  -0.86300 10.38200 21.00400
N  1.30900 9.79300 27.10700
C  1.22300 10.82400 26.23200
C  1.50800 12.06900 26.98600
C  1.50300 11.75700 28.31400
C  1.42300 10.34100 28.32600
C  1.55300 13.42200 26.29600
C  1.55800 12.26300 29.69200
O  1.68500 13.41500 30.08900
C  1.44800 11.03500 30.57500
C  2.56000 11.13900 31.53300
O  3.76600 10.89700 31.38900
O  2.04800 11.65700 32.70500
C  2.99200 11.90900 33.77200
H  2.69600 4.26200 29.63100
H  1.45200 4.99400 23.32300
H  0.89300 11.33000 24.16800
H  2.21800 8.40600 31.94400
H  1.06000 5.77700 31.61800
H  4.00400 5.85500 31.32000
H  2.97100 4.99800 32.29300
H  3.38200 6.77100 32.67900
H  -0.56100 7.28300 31.53300
H  -0.37400 8.97400 31.55200
H  -1.02500 7.59600 33.78000
H  -0.31800 9.19800 33.64500
H  3.22100 2.42800 28.75200
H  3.95300 1.82100 27.27400
H  2.28900 1.47200 27.52500
H  1.81900 0.77800 25.35900
H  3.66100 1.02700 25.13400
H  2.61600 0.48900 23.88200
H  -0.10600 7.10300 22.13500
H  1.34500 9.68200 22.18000
H  2.33200 8.19000 20.94400
H  1.31000 6.77100 20.73200
H  2.77000 6.84600 21.88100
H  -1.17200 10.51700 22.99600
H  -1.64500 8.88600 22.34300
H  -1.71200 9.95300 20.47200
H  -0.02300 10.09700 20.37100
H  -0.94100 11.46900 21.03900
H  1.59500 13.38900 25.20700
H  2.38500 14.03800 26.63500
H  0.58900 13.91000 26.44400
H  0.47700 11.06900 31.06800
H  3.81000 12.49700 33.35500
H  3.40600 11.00100 34.21100
H  2.51900 12.68100 34.37800
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


