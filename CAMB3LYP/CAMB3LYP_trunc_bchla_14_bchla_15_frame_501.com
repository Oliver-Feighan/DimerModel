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


