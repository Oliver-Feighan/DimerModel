%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 45.00000 3.44200 47.03300
C  42.74800 5.98100 46.47100
C  42.51500 1.08600 46.99500
C  47.33500 0.85800 47.41000
C  47.61400 5.55900 46.25100
N  42.89400 3.58700 46.53800
C  42.09600 4.73400 46.51200
C  40.62300 4.38600 46.28300
C  40.66000 2.87900 46.64300
C  42.13500 2.42800 46.65300
C  39.85400 2.54000 47.91600
C  40.17300 4.67400 44.83700
C  38.76200 5.33200 44.65300
H  38.03300 4.97800 43.41600
N  44.91300 1.33100 47.23000
C  43.80200 0.56100 47.28100
C  44.07100 -0.78100 47.70300
C  45.44900 -0.85000 47.96600
C  45.98000 0.51500 47.55400
C  42.99500 -1.83600 47.96800
C  46.25300 -2.09300 48.49400
O  45.64300 -3.18400 48.63300
C  47.77600 -2.03300 48.94400
N  47.15200 3.17200 46.56100
C  47.86200 2.08300 46.94400
C  49.35200 2.30000 46.87800
C  49.47200 3.78200 46.46200
C  48.00700 4.22600 46.45700
C  50.12800 1.90800 48.17900
C  50.29800 3.99400 45.10500
C  51.71400 4.54900 45.29800
N  45.19600 5.38500 46.52800
C  46.33300 6.17700 46.25500
C  45.93600 7.59000 46.19700
C  44.55800 7.53700 46.35400
C  44.18300 6.18500 46.42000
C  46.76000 8.81000 45.98900
C  43.34800 8.37000 46.40400
O  43.19200 9.59500 46.53700
C  42.11400 7.38900 46.52800
C  41.34500 7.62500 47.82100
O  41.80700 7.77900 48.95600
O  40.04700 7.69400 47.41000
C  39.07000 8.09100 48.41000
H  41.73100 0.32600 46.99100
H  48.15700 0.16300 47.59400
H  48.39400 6.28300 46.00700
H  39.99900 4.98400 46.94800
H  40.15100 2.41000 45.80100
H  39.07600 3.22300 48.25700
H  40.58100 2.53900 48.72800
H  39.36600 1.56800 47.84800
H  40.24600 3.80100 44.18900
H  40.97000 5.28600 44.41400
H  38.87100 6.41600 44.68600
H  38.11500 4.99000 45.46200
H  43.14300 -2.15000 49.00100
H  42.99500 -2.57100 47.16300
H  42.06000 -1.27600 47.95400
H  48.36300 -2.15900 48.03500
H  47.92900 -2.88500 49.60600
H  47.87200 -1.14900 49.57500
H  49.74800 1.66300 46.08700
H  49.86100 4.46700 47.21500
H  50.84800 2.67500 48.46300
H  50.67000 0.96800 48.07700
H  49.41000 1.70300 48.97200
H  49.82900 4.72000 44.44100
H  50.30900 3.07000 44.52700
H  52.39600 4.06600 44.59800
H  52.05000 4.35300 46.31600
H  51.79100 5.62600 45.15000
H  47.80500 8.50000 45.98500
H  46.42200 9.60500 46.65400
H  46.58100 9.11700 44.95800
H  41.44200 7.45600 45.67200
H  39.19800 7.59700 49.37400
H  38.01800 7.98300 48.14600
H  39.14800 9.16300 48.59300
Mg 6.73600 57.36100 41.87100
C  5.57600 54.03200 41.84400
C  9.54700 56.26000 40.29200
C  7.50300 60.62700 41.21800
C  3.48900 58.25400 43.03700
N  7.42000 55.35200 41.01200
C  6.89900 54.10900 41.23100
C  7.68800 53.00900 40.50000
C  9.08300 53.80000 40.39600
C  8.66300 55.23700 40.58100
C  10.16800 53.39100 41.42800
C  7.12500 52.63700 39.05500
C  7.43400 51.26500 38.47600
H  8.38100 51.07500 37.27200
N  8.34400 58.41300 40.90600
C  9.44600 57.74300 40.46000
C  10.47100 58.62900 40.12400
C  9.93400 59.94200 40.40300
C  8.52300 59.74000 40.83800
C  11.84700 58.22300 39.61000
C  10.74100 61.20700 40.21700
O  11.93000 61.13100 39.81300
C  10.22200 62.63600 40.60600
N  5.63100 59.10700 42.11600
C  6.19700 60.37100 41.79700
C  5.36400 61.47100 42.39400
C  3.99000 60.81500 42.65100
C  4.37000 59.30600 42.68800
C  5.80500 62.13900 43.71600
C  2.94100 61.09300 41.45000
C  2.15100 62.39300 41.46200
N  4.97100 56.31000 42.51500
C  3.75500 56.80200 42.96500
C  2.84300 55.78700 43.18500
C  3.53400 54.62400 42.80400
C  4.81800 55.00700 42.38100
C  1.45400 55.97300 43.70900
C  3.46600 53.25900 42.78000
O  2.61500 52.50900 43.21600
C  4.68000 52.84700 41.92100
C  4.09400 52.40300 40.61400
O  3.41000 53.09000 39.92100
O  4.40200 51.06200 40.37800
C  3.73000 50.54700 39.23600
H  10.54100 55.93300 39.97900
H  7.77500 61.68500 41.22200
H  2.57200 58.59300 43.52300
H  7.89200 52.14900 41.13700
H  9.49500 53.56500 39.41500
H  10.79300 52.70800 40.85200
H  9.75000 52.84400 42.27300
H  10.76800 54.24600 41.74000
H  7.52400 53.32000 38.30500
H  6.05400 52.83900 39.08200
H  6.45300 50.86600 38.21400
H  7.81100 50.65000 39.29300
H  12.72600 58.79500 39.90900
H  11.62200 58.24700 38.54400
H  12.11100 57.23800 39.99500
H  10.03600 62.67400 41.67900
H  9.32800 62.73100 39.99000
H  10.83200 63.47600 40.27400
H  5.23800 62.21700 41.61000
H  3.54200 61.14600 43.58800
H  5.15700 61.94700 44.57100
H  6.23700 63.12400 43.53800
H  6.70300 61.58800 43.99500
H  2.28600 60.23100 41.57500
H  3.43700 61.03400 40.48100
H  2.26900 62.89100 40.50000
H  2.53500 63.03900 42.25200
H  1.08300 62.20800 41.58000
H  0.80400 55.41700 43.03400
H  1.12000 57.01000 43.65700
H  1.32300 55.57400 44.71500
H  5.17600 52.01400 42.41800
H  2.78900 50.09000 39.54100
H  4.23500 49.69300 38.78300
H  3.52200 51.31600 38.49400


