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


