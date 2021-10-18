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
Mg -9.26600 18.27700 42.96500
C  -5.81100 17.58000 42.98100
C  -8.47400 21.42900 42.13900
C  -12.47500 18.67000 42.48700
C  -9.75200 14.74100 43.14100
N  -7.34400 19.34600 42.56200
C  -6.06200 18.93000 42.69300
C  -5.04000 20.02200 42.45000
C  -5.98500 21.29000 42.45900
C  -7.37700 20.70000 42.36700
C  -5.82600 22.23200 43.72300
C  -4.08200 20.00700 41.18000
C  -4.36600 19.09100 39.96000
H  -3.96500 19.61900 38.53200
N  -10.31800 19.93500 42.50700
C  -9.83700 21.10100 42.09700
C  -10.97100 21.97200 41.68400
C  -12.19800 21.17200 41.95900
C  -11.67600 19.85500 42.32500
C  -10.82500 23.41900 41.29700
C  -13.61600 21.63300 41.74500
O  -13.81100 22.79100 41.38400
C  -14.89900 20.73300 41.88300
N  -10.94800 16.84400 42.71500
C  -12.12100 17.37600 42.70200
C  -13.18200 16.19300 42.82000
C  -12.25700 14.97500 42.77500
C  -10.85500 15.54600 42.86600
C  -14.16200 16.18100 44.00500
C  -12.30000 14.18000 41.36900
C  -12.88700 12.75700 41.40300
N  -8.10600 16.49700 42.98200
C  -8.41600 15.17700 43.17700
C  -7.19300 14.39200 43.48700
C  -6.09500 15.27800 43.29700
C  -6.75800 16.56600 43.10100
C  -7.09800 12.86900 43.72900
C  -4.68100 15.45100 43.26700
O  -3.71100 14.67400 43.28400
C  -4.42500 16.97400 43.13600
C  -3.61000 17.39300 44.36300
O  -2.45500 17.69200 44.51200
O  -4.50700 17.60300 45.41600
C  -3.95600 18.13500 46.70100
H  -8.18700 22.44600 41.86300
H  -13.54300 18.86100 42.35700
H  -9.93500 13.71100 43.45300
H  -4.44700 20.04900 43.36400
H  -5.75600 21.91400 41.59500
H  -5.11000 21.81600 44.43100
H  -6.79600 22.11100 44.20400
H  -5.58800 23.24500 43.39800
H  -3.08000 19.78600 41.54800
H  -3.99700 21.02300 40.79300
H  -5.44200 18.91700 39.97000
H  -3.77000 18.17900 39.94600
H  -11.77100 23.96000 41.28900
H  -10.25700 23.43000 40.36600
H  -10.11700 23.80900 42.02900
H  -15.08400 20.15000 42.78600
H  -14.85100 19.96100 41.11500
H  -15.81700 21.30800 41.76500
H  -13.87700 16.19700 41.98100
H  -12.39400 14.18200 43.51100
H  -14.02600 17.04700 44.65300
H  -14.09900 15.25100 44.57100
H  -15.20700 16.25900 43.70500
H  -11.30400 14.05300 40.94700
H  -12.88100 14.81200 40.69700
H  -12.44900 12.05100 40.69700
H  -13.92900 12.87300 41.10500
H  -12.71100 12.44000 42.43100
H  -7.98800 12.30200 43.45800
H  -6.94100 12.78400 44.80400
H  -6.18900 12.50600 43.25000
H  -3.84800 17.16600 42.23200
H  -2.93200 17.88100 46.97400
H  -4.58000 17.65400 47.45400
H  -4.05900 19.20000 46.91000


