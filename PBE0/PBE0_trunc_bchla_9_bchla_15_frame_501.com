%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 35.26200 1.75200 29.96200
C  32.85500 2.58300 32.32000
C  37.64200 1.80100 32.35100
C  37.51600 1.12500 27.56800
C  32.74800 1.88400 27.52800
N  35.16600 2.13300 32.11400
C  34.11900 2.35200 32.87300
C  34.63200 2.45300 34.38300
C  36.03300 2.95200 34.09900
C  36.35200 2.21400 32.76000
C  36.05200 4.49500 34.04900
C  34.37500 1.08800 35.15500
C  34.37600 1.20900 36.68200
H  34.84400 2.50900 37.31300
N  37.24000 1.32700 29.98300
C  38.04600 1.32500 31.11900
C  39.36800 0.74800 30.73100
C  39.37200 0.59900 29.27900
C  38.02100 1.10700 28.92100
C  40.46600 0.47300 31.73100
C  40.47100 0.00700 28.28100
O  40.29400 -0.12300 27.06500
C  41.78300 -0.42500 28.90200
N  35.17700 1.65000 27.89700
C  36.22800 1.39600 27.12200
C  35.89000 1.47100 25.60100
C  34.35600 1.46600 25.66700
C  34.03700 1.61900 27.12100
C  36.49900 2.69400 24.86000
C  33.56700 0.32100 24.88900
C  33.00100 -0.79800 25.70000
N  33.26100 2.34500 29.84600
C  32.36300 2.28700 28.83400
C  31.06300 2.71100 29.21500
C  31.19600 2.84400 30.63300
C  32.52400 2.60100 30.90000
C  29.78600 2.98300 28.46300
C  30.54600 3.13100 31.87700
O  29.40100 3.42500 32.14500
C  31.52900 2.83000 33.05300
C  31.55800 3.97100 33.97300
O  31.92600 5.08300 33.62600
O  30.95500 3.71000 35.23000
C  30.85200 4.85600 36.14700
H  38.39100 1.89200 33.13900
H  38.28200 0.96200 26.80800
H  32.00200 1.96400 26.73500
H  34.02600 3.17400 34.93000
H  36.78800 2.59500 34.79900
H  36.64600 4.74400 34.92900
H  35.05700 4.93600 34.12200
H  36.63900 4.75700 33.16900
H  35.14200 0.35900 34.89200
H  33.36400 0.79400 34.87200
H  34.98900 0.39000 37.05700
H  33.35800 1.00500 37.01300
H  40.05900 0.17200 32.69700
H  41.18400 1.27200 31.91300
H  41.03100 -0.39400 31.38900
H  42.08900 0.44800 29.47800
H  42.58200 -0.61900 28.18700
H  41.63400 -1.33400 29.48500
H  36.25600 0.59300 25.06900
H  34.05000 2.44500 25.30000
H  35.69900 3.06000 24.21600
H  37.40400 2.52000 24.28000
H  36.80800 3.41800 25.61400
H  34.19600 -0.03500 24.07300
H  32.66800 0.75100 24.44800
H  31.91500 -0.71000 25.66100
H  33.14600 -0.73200 26.77900
H  33.24300 -1.79400 25.32900
H  29.91400 2.90200 27.38300
H  29.37400 3.93500 28.79900
H  28.98400 2.25500 28.58700
H  31.23100 1.94400 33.61300
H  30.55000 4.60700 37.16400
H  30.10200 5.48600 35.66800
H  31.81400 5.36300 36.08200
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


