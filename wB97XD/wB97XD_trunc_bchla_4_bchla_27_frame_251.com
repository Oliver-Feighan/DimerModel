%nproc=24
%mem=175gb
#p wB97XD/Def2SVP td=(nstates=5)

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
Mg -5.36100 24.76000 27.02500
C  -3.64200 26.63800 29.43600
C  -6.40000 22.65900 29.54700
C  -6.94500 22.99100 24.73300
C  -4.07500 26.91800 24.56900
N  -5.09300 24.76200 29.22400
C  -4.33400 25.60400 30.03700
C  -4.64200 25.43000 31.57600
C  -5.37000 24.10100 31.51700
C  -5.66000 23.80300 29.98800
C  -4.43700 22.92100 32.05400
C  -5.65800 26.61300 31.90200
C  -5.45100 27.38800 33.27700
H  -4.92400 26.43700 34.32100
N  -6.35800 23.02700 27.11800
C  -6.68000 22.31300 28.19300
C  -7.48400 21.16300 27.87400
C  -7.81800 21.31800 26.51400
C  -6.97100 22.47200 26.02900
C  -7.92400 20.11100 28.82500
C  -8.85200 20.60400 25.66800
O  -9.21200 21.01500 24.55900
C  -9.41800 19.25800 26.12900
N  -5.43800 24.92500 24.91100
C  -6.26300 24.05500 24.23800
C  -6.35200 24.50400 22.78100
C  -5.50900 25.87000 22.72900
C  -5.01000 25.97700 24.14600
C  -5.90500 23.33400 21.83500
C  -6.23300 27.10700 22.22100
C  -5.47900 27.85000 21.07400
N  -4.10000 26.41500 26.99800
C  -3.63700 27.19900 25.89500
C  -2.75500 28.21500 26.34700
C  -2.74300 28.01700 27.73800
C  -3.53100 26.89000 28.05100
C  -2.03600 29.28900 25.59600
C  -2.22900 28.50400 28.97800
O  -1.46700 29.40000 29.25600
C  -2.89200 27.68300 30.12000
C  -1.83900 27.24500 31.07200
O  -0.96900 26.39800 30.89000
O  -2.04600 27.91300 32.21100
C  -1.18400 27.58800 33.32800
H  -6.90800 22.00000 30.25400
H  -7.48200 22.35700 24.02300
H  -3.70300 27.64900 23.84800
H  -3.68500 25.58900 32.07200
H  -6.31300 24.16400 32.05900
H  -3.44500 23.27200 32.33900
H  -4.33700 22.17800 31.26300
H  -4.99100 22.51200 32.89900
H  -6.69800 26.29000 31.93800
H  -5.59700 27.30300 31.06000
H  -6.41300 27.65600 33.71400
H  -4.98200 28.36000 33.12200
H  -9.00400 20.22200 28.93100
H  -7.47900 20.25200 29.80900
H  -7.70000 19.10000 28.48600
H  -9.79600 18.88600 25.17700
H  -10.18200 19.34000 26.90300
H  -8.70400 18.49400 26.43600
H  -7.35300 24.74000 22.42000
H  -4.63900 25.71000 22.09200
H  -6.65200 23.29700 21.04200
H  -5.78500 22.37600 22.34100
H  -5.04100 23.62000 21.23400
H  -6.46800 27.83800 22.99400
H  -7.10900 26.70500 21.71000
H  -4.60800 27.27900 20.75200
H  -5.00600 28.76500 21.43000
H  -6.13000 28.01900 20.21600
H  -2.85300 29.97800 25.38000
H  -1.55300 28.80400 24.74800
H  -1.32200 29.82700 26.21800
H  -3.47700 28.43300 30.65400
H  -0.13900 27.68100 33.03100
H  -1.23800 26.56800 33.71000
H  -1.41700 28.31200 34.10900


