%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
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


