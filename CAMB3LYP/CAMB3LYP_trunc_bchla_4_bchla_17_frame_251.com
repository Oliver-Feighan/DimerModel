%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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
Mg 29.10200 58.37000 41.70600
C  26.21100 56.85000 40.52500
C  31.00300 55.92400 39.98500
C  31.78200 60.36200 41.82300
C  27.09800 61.07100 42.76300
N  28.62400 56.74200 40.07100
C  27.36300 56.19900 40.01300
C  27.34200 54.94600 39.19500
C  28.87800 54.50900 39.24100
C  29.59900 55.79900 39.77000
C  29.21600 53.21300 40.13700
C  26.66600 55.21800 37.75400
C  26.06800 54.05700 37.13700
H  26.60500 53.76400 35.72000
N  31.13700 58.12700 41.18100
C  31.72600 57.03700 40.48100
C  33.15300 57.29600 40.39000
C  33.42500 58.57400 40.93600
C  32.09000 59.05300 41.42900
C  34.07200 56.38900 39.66000
C  34.85200 59.26600 40.92400
O  35.77500 58.62800 40.42000
C  35.15900 60.50100 41.76400
N  29.32700 60.43000 42.12600
C  30.57900 60.94600 42.21600
C  30.58100 62.56100 42.56700
C  29.07200 62.72600 43.08400
C  28.42300 61.36300 42.68500
C  31.60200 63.02200 43.67500
C  28.34100 63.90500 42.43900
C  27.75100 65.04300 43.30600
N  27.08200 58.74100 41.92000
C  26.43600 59.84700 42.39100
C  25.02600 59.57700 42.45300
C  24.91200 58.38300 41.64600
C  26.21500 57.96500 41.32600
C  23.95900 60.50800 42.95200
C  23.95200 57.50800 41.00600
O  22.77000 57.42900 41.15100
C  24.75000 56.37400 40.25200
C  24.55300 54.95400 40.82700
O  24.66200 54.66400 41.99900
O  24.13200 54.08300 39.80100
C  23.73600 52.76500 40.19600
H  31.53000 55.10300 39.49600
H  32.55200 61.12000 41.66600
H  26.51800 61.92400 43.12000
H  26.89400 54.10800 39.73000
H  29.20300 54.31600 38.21800
H  29.65600 52.49600 39.44400
H  28.33800 52.79800 40.63100
H  29.96300 53.53300 40.86400
H  27.50600 55.54800 37.14300
H  26.02700 56.07100 37.98400
H  24.99800 54.25700 37.09500
H  26.28800 53.16600 37.72500
H  33.60300 55.42900 39.44500
H  34.92200 56.13300 40.29300
H  34.37800 56.81900 38.70600
H  34.62000 61.39400 41.44600
H  36.23600 60.66600 41.71900
H  34.91600 60.25600 42.79800
H  30.70200 63.02600 41.58900
H  28.95600 62.73400 44.16800
H  32.36800 63.65500 43.22800
H  32.08100 62.14600 44.11400
H  31.02600 63.57600 44.41700
H  27.57000 63.56900 41.74600
H  28.92000 64.44500 41.69000
H  26.66500 65.11700 43.36900
H  28.07200 65.96200 42.81500
H  28.11700 64.98900 44.33200
H  23.07500 59.90100 43.14900
H  23.49000 61.12600 42.18700
H  24.26100 61.15900 43.77200
H  24.32100 56.40300 39.25000
H  23.95500 52.21700 39.27900
H  22.65500 52.67800 40.30400
H  24.24700 52.46200 41.10900


