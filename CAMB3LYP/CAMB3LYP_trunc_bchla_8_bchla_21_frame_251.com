%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg 15.77600 52.08200 25.40200
C  17.41000 50.69800 28.22000
C  13.20100 52.90800 27.42100
C  14.25400 53.22200 22.73700
C  18.25300 50.46100 23.37600
N  15.27500 51.74600 27.60200
C  16.22500 51.39100 28.54500
C  15.65700 51.73300 29.88000
C  14.15600 51.89200 29.63600
C  14.16000 52.23400 28.12500
C  13.39700 50.53700 29.89600
C  16.28100 53.03000 30.62900
C  16.58500 52.99000 32.15800
H  15.44800 52.47100 33.08500
N  13.93500 52.94500 25.15100
C  13.06900 53.20000 26.08700
C  11.87700 53.79900 25.46200
C  12.22000 54.11300 24.11200
C  13.51400 53.39600 23.93000
C  10.48000 54.11300 26.11100
C  11.50000 54.93700 23.04500
O  11.87100 55.02600 21.85900
C  10.19300 55.49000 23.36500
N  16.14700 51.69600 23.28800
C  15.39400 52.41700 22.36500
C  16.06600 52.28500 20.94400
C  17.38400 51.52800 21.24300
C  17.29200 51.23400 22.74200
C  15.13200 51.68100 19.89500
C  18.69000 52.21600 20.77200
C  19.80100 51.32100 20.17500
N  17.47400 50.82500 25.67600
C  18.34200 50.20900 24.78000
C  19.41200 49.59700 25.52600
C  19.05500 49.71600 26.87000
C  17.86300 50.48800 26.94600
C  20.58200 49.01900 25.00900
C  19.46500 49.44600 28.25600
O  20.44600 48.90400 28.74000
C  18.40300 50.10900 29.19300
C  17.85700 49.01900 30.09800
O  16.99100 48.18400 29.96600
O  18.53900 49.25300 31.25700
C  18.03900 48.60700 32.49400
H  12.40000 53.26000 28.07500
H  13.81800 53.79200 21.91400
H  19.12700 50.15700 22.79700
H  15.77300 50.97700 30.65700
H  13.74100 52.77300 30.12500
H  12.38500 50.66400 30.27900
H  13.84800 49.93600 30.68600
H  13.32800 50.06400 28.91700
H  15.69100 53.91000 30.37500
H  17.23100 53.29600 30.16400
H  16.68900 54.04900 32.39700
H  17.49400 52.38900 32.18400
H  10.24800 55.17100 26.23100
H  10.47100 53.63100 27.08900
H  9.61500 53.80800 25.52200
H  10.23700 56.09500 24.27100
H  9.45700 54.69600 23.49200
H  9.96700 56.22700 22.59500
H  16.38500 53.29100 20.67100
H  17.33800 50.51100 20.85400
H  15.24100 50.59600 19.90400
H  15.36700 51.95500 18.86700
H  14.10200 51.79300 20.23200
H  19.14000 52.87200 21.51700
H  18.48100 52.96500 20.00700
H  20.48500 51.15800 21.00800
H  20.30400 51.95600 19.44500
H  19.37300 50.39500 19.79200
H  21.13200 48.28300 25.59500
H  21.23200 49.84300 24.71300
H  20.26200 48.46700 24.12600
H  18.86400 50.90900 29.77200
H  18.09100 49.38900 33.25100
H  18.61600 47.75900 32.86200
H  16.99700 48.31500 32.36600


