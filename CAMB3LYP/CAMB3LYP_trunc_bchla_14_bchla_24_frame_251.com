%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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


