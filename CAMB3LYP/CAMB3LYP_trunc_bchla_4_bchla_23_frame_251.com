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
Mg -9.42700 40.51200 41.96000
C  -8.27600 37.27100 41.01300
C  -7.18900 41.70700 39.82100
C  -11.06300 43.46900 41.99400
C  -12.24300 38.82400 43.34600
N  -7.88600 39.63200 40.71800
C  -7.50800 38.28500 40.52900
C  -6.30600 38.11700 39.54700
C  -5.96100 39.60100 39.23600
C  -7.01600 40.39000 40.04000
C  -4.48700 39.82100 39.61500
C  -6.54600 37.26600 38.22600
C  -5.77300 35.95800 38.10300
H  -4.88600 35.77200 36.85500
N  -9.14300 42.33400 41.14400
C  -8.11800 42.62800 40.24500
C  -8.09400 44.00500 39.95700
C  -9.10700 44.61100 40.67700
C  -9.83500 43.52300 41.29800
C  -7.10100 44.80200 39.03200
C  -9.33100 46.10500 40.68600
O  -8.59600 46.89800 40.17700
C  -10.57800 46.61100 41.56500
N  -11.40100 41.07400 42.51100
C  -11.80600 42.38000 42.46200
C  -13.34000 42.40200 42.85400
C  -13.69600 41.00800 43.39700
C  -12.37200 40.22400 43.05500
C  -13.75900 43.56700 43.70500
C  -14.84900 40.42700 42.58800
C  -15.68100 39.35600 43.24200
N  -10.03000 38.45700 42.37900
C  -11.16900 37.95500 42.96200
C  -11.08800 36.54300 43.06100
C  -9.91600 36.15600 42.36900
C  -9.38800 37.38300 41.87700
C  -12.09900 35.62000 43.61300
C  -9.03800 35.08300 41.88300
O  -9.11400 33.86600 42.09900
C  -7.91300 35.74600 40.99900
C  -6.65000 35.63000 41.67800
O  -6.37600 36.35300 42.66500
O  -5.96400 34.63500 41.08200
C  -4.53900 34.42700 41.45900
H  -6.43900 42.13800 39.15400
H  -11.53000 44.45500 42.04300
H  -13.01900 38.25900 43.86600
H  -5.46700 37.65700 40.06900
H  -6.10300 39.87100 38.19000
H  -4.30300 40.89300 39.55000
H  -3.75300 39.31500 38.98800
H  -4.26800 39.50000 40.63400
H  -6.11700 37.88500 37.43800
H  -7.61400 37.12100 38.06300
H  -6.50000 35.14700 38.15800
H  -5.08700 35.86800 38.94600
H  -7.73100 45.35500 38.33500
H  -6.44100 44.13800 38.47400
H  -6.37700 45.43500 39.54500
H  -10.52000 47.69900 41.55500
H  -10.39300 46.15100 42.53600
H  -11.42100 46.18000 41.02500
H  -13.83700 42.54600 41.89400
H  -13.95800 40.94300 44.45300
H  -14.59900 43.25300 44.32400
H  -14.03900 44.35800 43.00900
H  -12.95300 44.01500 44.28500
H  -14.53800 40.03800 41.61800
H  -15.54600 41.19600 42.25500
H  -15.47900 39.25300 44.30800
H  -15.45200 38.41100 42.74900
H  -16.73600 39.48600 43.00100
H  -12.69800 36.09600 44.39000
H  -11.53400 34.78400 44.02600
H  -12.74800 35.36800 42.77400
H  -7.89000 35.41100 39.96200
H  -4.47300 33.34300 41.35500
H  -4.28400 34.79600 42.45200
H  -3.86400 34.76500 40.67200


