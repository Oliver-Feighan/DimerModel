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
Mg -0.43700 44.18600 24.82500
C  1.43400 43.80100 27.72700
C  -3.18200 43.12800 26.59100
C  -2.25600 44.31000 21.89400
C  2.46600 44.70300 22.96600
N  -0.82100 43.63800 26.94100
C  0.04100 43.65300 27.97200
C  -0.60100 43.41700 29.35900
C  -2.06400 42.96100 28.94000
C  -2.04900 43.21000 27.35600
C  -2.28200 41.47000 29.38200
C  -0.53600 44.82000 30.15500
C  0.12400 44.79000 31.56500
H  0.16600 43.36300 32.22200
N  -2.40500 43.72200 24.25200
C  -3.41200 43.37700 25.17700
C  -4.69200 43.40700 24.50200
C  -4.45000 43.69100 23.15200
C  -2.97800 44.00900 23.04400
C  -5.95900 42.85900 25.17400
C  -5.47200 43.73000 21.98500
O  -5.19700 43.99600 20.80000
C  -6.95700 43.49800 22.31300
N  0.09900 44.39300 22.73100
C  -0.90800 44.46700 21.72100
C  -0.26100 44.51400 20.32100
C  1.23600 44.92600 20.70600
C  1.26600 44.70400 22.24100
C  -0.47400 43.18300 19.61900
C  1.52900 46.44600 20.35700
C  0.65700 47.66400 20.81600
N  1.53600 44.17200 25.19400
C  2.63200 44.38000 24.35100
C  3.83100 44.31600 25.12200
C  3.42700 44.13900 26.43300
C  2.02300 43.94900 26.43800
C  5.25700 44.36900 24.61300
C  3.89700 44.01900 27.75600
O  5.03000 43.96500 28.11500
C  2.62100 43.67600 28.63100
C  2.82000 42.26800 29.12000
O  2.22200 41.30700 28.68300
O  3.77300 42.30100 30.05900
C  4.37700 41.01200 30.29200
H  -4.07700 42.87100 27.16100
H  -2.79100 44.11000 20.96400
H  3.38600 44.85200 22.39600
H  -0.09900 42.67700 29.98200
H  -2.86900 43.57400 29.34600
H  -2.77400 40.83600 28.64400
H  -2.90900 41.64900 30.25500
H  -1.32000 41.10700 29.74300
H  -1.54200 45.23900 30.14700
H  0.04200 45.50300 29.53200
H  -0.52900 45.47000 32.11100
H  1.13900 45.18400 31.60700
H  -6.35300 42.07900 24.52200
H  -6.71900 43.64100 25.20200
H  -5.72100 42.49500 26.17300
H  -7.06800 42.43500 22.53100
H  -7.49100 43.66100 21.37700
H  -7.30300 44.26900 23.00100
H  -0.72000 45.37600 19.83600
H  1.94200 44.36800 20.09200
H  -1.38200 42.67100 19.93800
H  0.36900 42.51600 19.79800
H  -0.61300 43.32900 18.54800
H  1.65900 46.51600 19.27700
H  2.46900 46.62800 20.87800
H  -0.12600 47.82100 20.07400
H  1.34000 48.51400 20.83800
H  0.27600 47.62500 21.83700
H  5.47400 45.30000 24.08900
H  5.31100 43.53800 23.91000
H  6.03200 44.40300 25.38000
H  2.44000 44.27600 29.52300
H  5.42200 41.06200 30.59700
H  4.25200 40.18600 29.59200
H  3.93900 40.47400 31.13300


