%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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
Mg 46.35100 44.49200 43.30700
C  43.07600 43.62000 43.02100
C  47.22500 41.16500 42.48500
C  49.41300 45.58500 42.54300
C  45.37900 47.95800 43.59800
N  45.25500 42.66900 42.85200
C  43.97900 42.48600 42.76800
C  43.57300 40.97200 42.67800
C  44.93400 40.18700 42.80300
C  45.92900 41.39700 42.74700
C  45.03400 39.17700 43.94700
C  42.68400 40.64700 41.33600
C  42.79700 41.58100 40.14200
H  42.36600 41.11900 38.75000
N  48.07200 43.49800 42.60100
C  48.26200 42.15300 42.39300
C  49.56000 41.92800 41.92400
C  50.24000 43.20600 41.88200
C  49.22200 44.19700 42.32000
C  50.08500 40.51700 41.73200
C  51.70300 43.39300 41.44100
O  52.34200 42.37400 41.09500
C  52.40800 44.64100 41.38700
N  47.26200 46.56500 42.93800
C  48.55900 46.63600 42.83800
C  49.10000 48.08900 42.96000
C  47.81700 48.85400 43.35200
C  46.67900 47.75500 43.33200
C  50.17200 48.26600 43.97400
C  47.52200 50.10400 42.39600
C  47.68100 51.53700 43.05200
N  44.63200 45.62300 43.43200
C  44.35200 46.98700 43.58200
C  42.86900 47.16200 43.67400
C  42.37400 45.87600 43.41500
C  43.46600 44.97500 43.24000
C  42.10000 48.37100 43.90000
C  41.16700 45.21000 43.23100
O  39.97800 45.55600 43.29100
C  41.54900 43.71200 43.06200
C  40.81500 42.90200 44.15700
O  40.09100 41.93300 44.00400
O  41.06300 43.39200 45.42100
C  40.39300 42.73400 46.54900
H  47.39200 40.13000 42.18000
H  50.42700 45.98200 42.46400
H  45.10600 48.98300 43.85700
H  43.05600 40.74200 43.61000
H  45.17500 39.63500 41.89400
H  45.78600 39.50800 44.66300
H  45.34200 38.19000 43.60000
H  44.10300 38.96300 44.47200
H  41.68100 40.58600 41.75900
H  42.91700 39.70400 40.84100
H  43.85800 41.82900 40.12800
H  42.26900 42.48700 40.44000
H  50.98600 40.19700 42.25600
H  50.18000 40.32500 40.66300
H  49.31900 39.81000 42.05000
H  52.54000 45.19200 42.31800
H  51.83700 45.37400 40.81800
H  53.35500 44.58600 40.85100
H  49.55600 48.55200 42.08500
H  47.80800 49.17400 44.39400
H  51.17700 48.38000 43.56800
H  50.31400 47.34900 44.54700
H  49.94700 49.11900 44.61400
H  46.47000 50.06600 42.11000
H  48.18600 50.18000 41.53600
H  48.24700 52.18400 42.38300
H  48.32000 51.43300 43.92900
H  46.66400 51.92200 43.12900
H  42.74600 48.98700 44.52600
H  41.17900 48.14600 44.43900
H  41.83000 48.74800 42.91400
H  41.05300 43.23500 42.21700
H  40.32500 43.44800 47.37000
H  40.87800 41.87000 47.00300
H  39.41700 42.31700 46.29900


