%nproc=24
%mem=175gb
#p wB97XD/Def2SVP td=(nstates=5)

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


