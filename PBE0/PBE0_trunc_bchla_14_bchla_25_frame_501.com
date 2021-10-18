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
Mg -2.33100 34.27800 27.06700
C  -3.29300 32.73500 29.99200
C  -0.77900 36.82100 29.04000
C  -1.68900 35.99300 24.34500
C  -4.00500 31.86100 25.19100
N  -2.10600 34.64700 29.26400
C  -2.54800 33.87100 30.25800
C  -2.04400 34.45500 31.61800
C  -1.35700 35.84600 31.23800
C  -1.46000 35.80800 29.70100
C  -2.07400 37.11000 31.67200
C  -1.16200 33.50400 32.44100
C  -1.65900 33.38600 33.87200
H  -0.69900 33.03500 35.01800
N  -1.29600 36.10900 26.75300
C  -0.70200 36.96400 27.63600
C  -0.06500 38.02800 26.99800
C  -0.40300 37.90700 25.60500
C  -1.15900 36.61100 25.50300
C  0.74500 39.09500 27.72000
C  -0.12100 38.87200 24.46500
O  -0.65300 38.65500 23.35900
C  0.65800 40.10800 24.66200
N  -2.76200 33.93100 25.11000
C  -2.33900 34.76200 24.09300
C  -2.78500 34.31300 22.65200
C  -3.27700 32.84000 22.98500
C  -3.37800 32.87200 24.52100
C  -3.86900 35.12000 21.83400
C  -2.50400 31.61000 22.32400
C  -1.18700 31.25300 22.96100
N  -3.52800 32.72900 27.43000
C  -4.16000 31.81100 26.61300
C  -4.85300 30.83200 27.34600
C  -4.49800 31.15000 28.66200
C  -3.68000 32.31700 28.68800
C  -5.57400 29.72900 26.69600
C  -4.80700 30.73800 29.97700
O  -5.59800 29.85900 30.32500
C  -4.00900 31.74000 30.96400
C  -4.95800 32.38400 31.85500
O  -5.73200 33.27000 31.57300
O  -4.95600 31.78800 33.04100
C  -5.91200 32.19500 34.02600
H  -0.35300 37.61900 29.65100
H  -1.55300 36.44800 23.36200
H  -4.51800 31.12700 24.56600
H  -2.99900 34.65700 32.10200
H  -0.41300 35.85600 31.78400
H  -2.94800 36.94200 32.30200
H  -2.45800 37.64200 30.80200
H  -1.35600 37.81800 32.08700
H  -0.10200 33.75800 32.45800
H  -1.18900 32.52900 31.95500
H  -2.34600 32.53900 33.87800
H  -2.14000 34.33700 34.10000
H  0.94800 38.81600 28.75400
H  0.13200 39.99100 27.82000
H  1.62000 39.38400 27.13900
H  0.44000 40.83300 23.87700
H  1.70900 39.83900 24.55200
H  0.43200 40.54800 25.63300
H  -1.93500 34.12600 21.99500
H  -4.26000 32.69700 22.53800
H  -4.17100 35.97400 22.43900
H  -4.79700 34.54700 21.82100
H  -3.62200 35.43900 20.82100
H  -2.25200 31.92600 21.31200
H  -3.11500 30.71700 22.19800
H  -1.06000 30.17800 22.83800
H  -0.96900 31.52500 23.99300
H  -0.42200 31.79500 22.40500
H  -5.20900 29.35700 25.73900
H  -6.66000 29.81200 26.67700
H  -5.43800 28.93700 27.43200
H  -3.30800 31.14200 31.54600
H  -5.54400 33.01100 34.64900
H  -6.05600 31.34600 34.69400
H  -6.83100 32.52000 33.53800


