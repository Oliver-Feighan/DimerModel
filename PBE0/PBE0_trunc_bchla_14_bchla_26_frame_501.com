%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg -9.66500 18.33400 42.78400
C  -6.24000 17.80500 42.49600
C  -9.09900 21.67600 42.05300
C  -12.91200 18.87400 42.48300
C  -10.14800 14.91400 43.00500
N  -7.91600 19.56100 42.17100
C  -6.59800 19.16200 42.25800
C  -5.63300 20.29600 42.07600
C  -6.57300 21.55800 42.31900
C  -7.95100 20.93600 42.13600
C  -6.37100 22.54100 43.56800
C  -4.81400 20.28500 40.70300
C  -5.51800 19.69700 39.44400
H  -5.27900 20.48100 38.17600
N  -10.85700 20.01500 42.21400
C  -10.42000 21.32200 41.96500
C  -11.51000 22.18000 41.58900
C  -12.66600 21.36900 41.68400
C  -12.20600 20.01400 42.11100
C  -11.29500 23.68900 41.46200
C  -14.07900 21.93000 41.40500
O  -14.25800 23.08100 40.99400
C  -15.20100 21.03000 41.71200
N  -11.24800 17.01300 42.63400
C  -12.52300 17.52900 42.79400
C  -13.59000 16.43700 43.02800
C  -12.64200 15.17100 42.94600
C  -11.28500 15.70200 42.80600
C  -14.41600 16.47100 44.27500
C  -12.97300 14.21900 41.68800
C  -13.02000 12.77100 42.10100
N  -8.45700 16.64500 42.95800
C  -8.80300 15.30900 43.04100
C  -7.58100 14.47400 43.08900
C  -6.53000 15.43300 42.87800
C  -7.12900 16.69600 42.80500
C  -7.42000 13.00900 43.03700
C  -5.08700 15.65900 42.64900
O  -4.19400 14.85200 42.53200
C  -4.84400 17.18800 42.43700
C  -3.86600 17.68600 43.40500
O  -2.75700 18.18900 43.22300
O  -4.28400 17.22800 44.68200
C  -3.45300 17.67600 45.82000
H  -8.86100 22.74200 42.07500
H  -13.99400 18.96100 42.60100
H  -10.40600 13.86400 43.15400
H  -4.85200 20.25200 42.83500
H  -6.35600 22.13700 41.42100
H  -5.49500 22.34800 44.18700
H  -7.18200 22.50000 44.29500
H  -6.14100 23.54200 43.20200
H  -3.84000 19.81700 40.84200
H  -4.55300 21.32400 40.50200
H  -6.59000 19.59400 39.61500
H  -5.14100 18.69400 39.24600
H  -11.94900 24.16700 40.73200
H  -10.29700 23.87500 41.06600
H  -11.49900 24.10000 42.45100
H  -15.27200 20.92300 42.79400
H  -15.05600 20.07500 41.20600
H  -16.08400 21.54300 41.33100
H  -14.18800 16.36800 42.11900
H  -12.69300 14.77700 43.96100
H  -14.92900 15.51000 44.31600
H  -15.20300 17.20800 44.11500
H  -13.79200 16.69500 45.14100
H  -12.21000 14.42700 40.93800
H  -13.86700 14.45400 41.10900
H  -12.15800 12.27100 41.66000
H  -13.81800 12.22200 41.60000
H  -12.89800 12.69100 43.18100
H  -8.34400 12.43000 43.00100
H  -6.93100 12.73100 43.97000
H  -6.96900 12.89900 42.05100
H  -4.46800 17.26400 41.41700
H  -4.21000 17.98200 46.54300
H  -2.74000 18.49300 45.70900
H  -3.06600 16.75200 46.25000


