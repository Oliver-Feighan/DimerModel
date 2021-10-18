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
Mg 0.26000 43.34200 24.75600
C  2.25300 43.09900 27.61800
C  -2.44200 42.42000 26.69400
C  -1.57400 42.87700 21.91500
C  3.14600 43.60500 22.86000
N  -0.01100 42.92700 26.88700
C  0.88400 43.10000 27.86800
C  0.29200 43.09400 29.28500
C  -1.20300 42.68500 28.99500
C  -1.29800 42.75900 27.37000
C  -1.68200 41.32300 29.54200
C  0.44600 44.41100 30.04600
C  0.46000 44.23200 31.57500
H  0.30000 42.82500 32.16800
N  -1.74500 42.84700 24.34000
C  -2.74500 42.55800 25.29700
C  -3.99500 42.34000 24.59700
C  -3.74800 42.40500 23.17500
C  -2.33100 42.75300 23.08200
C  -5.28700 42.18000 25.35000
C  -4.75800 42.32200 21.94400
O  -4.44000 42.53400 20.79900
C  -6.18400 41.93600 22.27200
N  0.72800 43.16000 22.71000
C  -0.18900 43.01100 21.71100
C  0.47200 43.04900 20.28400
C  1.92800 43.34600 20.63900
C  1.92700 43.36400 22.15600
C  0.26000 41.87700 19.24300
C  2.38800 44.69400 20.03600
C  2.04300 46.00200 20.78200
N  2.32100 43.50900 25.09100
C  3.38200 43.61100 24.22900
C  4.61600 43.69300 24.93600
C  4.22700 43.48400 26.33500
C  2.89200 43.36300 26.33900
C  5.99200 43.88800 24.47400
C  4.71200 43.25500 27.68900
O  5.85600 43.34100 28.14800
C  3.47500 43.00100 28.60000
C  3.56600 41.71300 29.28000
O  3.33400 40.61600 28.78300
O  4.17100 41.89000 30.49900
C  4.75700 40.73100 31.20000
H  -3.23100 42.10200 27.37900
H  -2.18100 42.91300 21.00800
H  4.04900 43.78500 22.27300
H  0.67700 42.26400 29.87700
H  -1.92300 43.45600 29.26700
H  -2.57300 41.63600 30.08600
H  -0.94100 40.89700 30.21800
H  -1.87200 40.62800 28.72400
H  -0.32500 45.13200 29.77700
H  1.38300 44.89600 29.77100
H  -0.37400 44.87500 31.85700
H  1.32100 44.72000 32.03400
H  -5.00100 42.06900 26.39600
H  -5.81000 41.30000 24.97500
H  -5.96200 43.02800 25.22500
H  -6.37400 40.88800 22.50600
H  -6.81500 42.27400 21.45000
H  -6.38700 42.55900 23.14300
H  -0.03500 43.91100 19.85100
H  2.51500 42.52100 20.23700
H  1.17400 41.51200 18.77600
H  -0.34600 42.26100 18.42300
H  -0.30900 41.09200 19.74200
H  1.97300 44.71300 19.02900
H  3.46900 44.59400 19.93400
H  1.75200 46.69500 19.99300
H  2.89700 46.44600 21.29300
H  1.29400 45.88900 21.56500
H  6.47500 44.58100 25.16200
H  6.09900 44.26000 23.45500
H  6.57000 42.97800 24.63600
H  3.48200 43.84300 29.29200
H  4.58400 39.79600 30.66600
H  4.26600 40.64900 32.17000
H  5.80700 41.02000 31.22900


