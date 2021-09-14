%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 1.27700 7.97600 26.81300
C  1.70500 10.25900 29.27500
C  1.82700 5.41400 29.09300
C  0.85500 5.69300 24.34900
C  1.21800 10.51800 24.46400
N  1.59300 7.92000 28.89400
C  1.77500 8.96200 29.74400
C  1.81600 8.52800 31.18600
C  2.13100 6.98700 31.08800
C  1.79700 6.74800 29.57900
C  3.53500 6.42400 31.56400
C  0.41400 8.75400 31.91200
C  0.55800 8.89600 33.45600
H  1.71800 8.24400 34.22700
N  1.40400 5.81900 26.67300
C  1.63700 4.96400 27.78200
C  1.59900 3.60000 27.34800
C  1.32800 3.64600 25.93400
C  1.17100 5.11200 25.62400
C  1.72200 2.42000 28.24700
C  1.18400 2.49900 24.94000
O  1.05900 2.62000 23.68600
C  1.20700 1.06700 25.44700
N  1.13600 8.07200 24.71300
C  0.77100 6.99400 23.90500
C  0.28500 7.38700 22.48000
C  0.38900 9.00800 22.61900
C  0.90800 9.24400 23.99000
C  1.23400 6.77600 21.44100
C  -0.88200 9.73900 22.30100
C  -0.72400 10.77600 21.22600
N  1.52400 9.91600 26.73700
C  1.45400 10.87900 25.73800
C  1.65900 12.16800 26.28700
C  1.81400 12.01500 27.69100
C  1.62000 10.64100 27.94600
C  1.67300 13.48200 25.60500
C  2.15100 12.69400 28.92300
O  2.45500 13.83200 29.22300
C  2.00700 11.52600 30.01000
C  3.06400 11.59000 31.05300
O  4.23800 11.79300 30.87600
O  2.59100 11.43200 32.29900
C  3.63300 11.47700 33.36400
H  2.18200 4.70500 29.84400
H  0.85400 4.95200 23.54700
H  1.09900 11.38200 23.80700
H  2.67100 8.94100 31.72100
H  1.35600 6.44600 31.63000
H  4.07900 7.25400 32.01500
H  4.16900 5.94900 30.81700
H  3.32200 5.63000 32.28000
H  -0.18400 7.89400 31.61000
H  -0.01700 9.66700 31.50100
H  -0.39200 8.49700 33.81200
H  0.49500 9.96300 33.66800
H  0.81200 2.12400 28.76900
H  2.52200 2.67200 28.94200
H  2.18600 1.56000 27.76500
H  2.26000 0.80200 25.54100
H  0.75700 0.48100 24.64600
H  0.71900 0.92900 26.41200
H  -0.69500 6.92500 22.36000
H  1.13100 9.36700 21.90500
H  2.13800 6.31200 21.83600
H  1.41900 7.51100 20.65800
H  0.69100 5.93800 21.00300
H  -1.27500 10.30600 23.14500
H  -1.71400 9.09900 22.00700
H  -1.48700 10.70900 20.45100
H  0.24000 10.68600 20.72600
H  -0.68400 11.71900 21.77200
H  2.02800 13.48300 24.57400
H  2.51200 13.99500 26.07600
H  0.69500 13.95600 25.68900
H  1.08300 11.88900 30.46100
H  3.95200 12.49700 33.58300
H  4.44500 10.77600 33.17000
H  3.18800 11.10700 34.28800
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


