%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 8.88300 3.35600 27.86700
C  10.38300 1.75300 30.56400
C  7.38700 5.53100 30.05200
C  7.19500 4.67800 25.21800
C  10.33700 0.83700 25.82400
N  8.85600 3.59100 30.04700
C  9.65000 2.88400 30.91100
C  9.43500 3.49700 32.29700
C  8.27700 4.46700 32.14400
C  8.19100 4.56300 30.62500
C  6.95700 4.16900 32.80300
C  10.74000 4.12300 32.98500
C  10.49500 4.94200 34.29300
H  10.56500 4.20300 35.65500
N  7.60700 4.98100 27.67200
C  7.08500 5.72600 28.67800
C  6.18300 6.66800 28.06800
C  6.13000 6.40800 26.66700
C  6.96700 5.31400 26.45500
C  5.47200 7.68700 28.90400
C  5.19000 7.09800 25.63600
O  4.95300 6.60000 24.55000
C  4.43400 8.30700 25.92000
N  8.80500 2.84800 25.93900
C  8.07600 3.55600 24.98400
C  8.26300 2.99400 23.56400
C  9.45300 1.99700 23.77400
C  9.56600 1.92700 25.28800
C  7.04500 2.56000 22.85900
C  10.87100 2.27700 23.02500
C  11.45600 1.24800 22.05300
N  10.10500 1.54900 28.11900
C  10.57400 0.66100 27.19200
C  11.35500 -0.36000 27.80500
C  11.38100 0.07200 29.18600
C  10.52100 1.19800 29.31900
C  12.01900 -1.46000 27.06800
C  11.98300 -0.20000 30.47100
O  12.81500 -0.99400 30.86400
C  11.31100 0.93300 31.44800
C  10.73600 0.32300 32.67400
O  9.77600 -0.41600 32.67100
O  11.39000 0.69300 33.78900
C  10.74800 0.34700 35.11500
H  6.95000 6.26000 30.73800
H  6.63700 5.11600 24.38800
H  10.70900 0.18100 25.03400
H  9.14100 2.60500 32.85000
H  8.55100 5.42600 32.58300
H  6.38200 3.89500 31.91900
H  6.56900 5.01600 33.36900
H  6.99000 3.26400 33.41000
H  11.10800 4.77800 32.19500
H  11.45800 3.33300 33.20200
H  9.50400 5.39400 34.30500
H  11.21400 5.76100 34.31200
H  5.72000 7.75000 29.96400
H  4.42300 7.39300 28.89100
H  5.49900 8.68800 28.47300
H  3.85400 8.38500 26.83900
H  3.77700 8.64800 25.12000
H  5.21700 9.06500 25.91700
H  8.64800 3.90400 23.10300
H  8.95400 1.07600 23.47300
H  6.77200 3.11000 21.95800
H  6.21600 2.63700 23.56300
H  7.08300 1.53300 22.49600
H  11.62300 2.42400 23.80000
H  10.95600 3.19100 22.43600
H  12.28200 0.65800 22.45100
H  11.66900 1.63000 21.05500
H  10.70100 0.47500 21.90900
H  11.30800 -1.87500 26.35400
H  12.43500 -2.21900 27.73000
H  12.73100 -1.00300 26.38000
H  12.14100 1.55900 31.77400
H  11.59600 0.40500 35.79700
H  10.32500 -0.65500 35.04600
H  10.03300 1.14100 35.33100
Mg 35.66000 1.74600 29.23200
C  33.47200 2.34100 31.89500
C  38.12100 1.24100 31.62500
C  37.79000 1.61500 26.77900
C  33.01700 2.53400 27.01000
N  35.87600 2.18500 31.54100
C  34.82500 2.24700 32.36500
C  35.34000 2.18200 33.76400
C  36.88000 2.27300 33.63100
C  37.00600 1.84600 32.18200
C  37.59500 3.59000 33.93800
C  34.89600 0.87200 34.48900
C  34.26200 0.97200 35.91400
H  34.68800 2.14100 36.85600
N  37.69800 1.31100 29.19900
C  38.45500 0.95400 30.30700
C  39.74800 0.60100 29.81800
C  39.73600 0.82500 28.38300
C  38.39000 1.25600 28.03600
C  40.94400 0.20000 30.70900
C  40.89400 0.72300 27.34900
O  40.78100 0.94900 26.18600
C  42.22500 0.23800 27.88900
N  35.39600 2.02600 27.20700
C  36.47200 1.98000 26.39000
C  36.13800 2.22700 24.89300
C  34.55100 2.25100 24.96900
C  34.26500 2.24500 26.46600
C  36.85300 3.59900 24.36100
C  33.86400 1.14100 24.11300
C  33.41700 -0.12600 24.82000
N  33.60700 2.29400 29.34300
C  32.66300 2.48800 28.37100
C  31.34700 2.76700 29.04600
C  31.64200 2.71200 30.42700
C  33.04500 2.46100 30.58700
C  30.05700 3.20300 28.40600
C  31.07900 2.83900 31.76900
O  29.91300 3.09600 32.08200
C  32.20100 2.44500 32.72500
C  32.26400 3.54000 33.75100
O  32.77700 4.63100 33.58200
O  31.49900 3.13900 34.75300
C  31.03100 4.13600 35.70300
H  38.92600 1.22300 32.36200
H  38.44500 1.39800 25.93300
H  32.19900 2.76200 26.32300
H  34.93300 3.02400 34.32300
H  37.28500 1.57100 34.36000
H  36.96900 4.41100 34.28800
H  38.13600 3.93900 33.05800
H  38.32400 3.36100 34.71500
H  35.66700 0.10200 34.50700
H  34.14500 0.53000 33.77700
H  34.54900 0.03000 36.38300
H  33.17800 1.00900 35.80900
H  41.84000 0.79100 30.51900
H  41.29800 -0.80900 30.49700
H  40.73900 0.22800 31.78000
H  42.13800 -0.68600 28.46100
H  42.66900 1.05000 28.46500
H  42.94900 0.14100 27.08000
H  36.48000 1.41700 24.24900
H  34.11700 3.15500 24.54200
H  36.02900 4.30300 24.24400
H  37.25500 3.51600 23.35100
H  37.53000 4.05700 25.08300
H  34.65800 0.79900 23.44900
H  33.05000 1.52900 23.50100
H  34.25200 -0.36700 25.47900
H  33.30700 -0.80900 23.97800
H  32.56000 -0.11800 25.49300
H  30.04900 3.41600 27.33700
H  29.71100 4.11800 28.88700
H  29.33300 2.43700 28.68200
H  31.98900 1.57300 33.34400
H  31.22900 5.19300 35.52800
H  31.60100 3.89200 36.60000
H  30.07700 3.79700 36.10700


