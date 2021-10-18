%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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


