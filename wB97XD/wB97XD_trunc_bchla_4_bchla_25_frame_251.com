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
Mg -2.51200 34.32000 26.71000
C  -3.70600 32.37000 29.49400
C  -1.09600 36.45200 28.85600
C  -2.23500 36.47200 24.15600
C  -4.45400 32.24200 24.65700
N  -2.47200 34.41100 28.88300
C  -2.84200 33.43600 29.82500
C  -2.27900 33.77800 31.25700
C  -1.51600 35.11400 31.00700
C  -1.64900 35.32300 29.47500
C  -1.96800 36.32500 31.78300
C  -1.39500 32.57400 31.88800
C  -1.81300 32.12300 33.29300
H  -0.66900 32.23700 34.36800
N  -1.63700 36.15800 26.50500
C  -1.11000 36.89500 27.48800
C  -0.61000 38.17800 26.99100
C  -1.04200 38.20400 25.62000
C  -1.69000 36.91400 25.36100
C  0.10900 39.14600 27.82900
C  -0.90500 39.45100 24.69100
O  -1.22300 39.41800 23.51400
C  -0.35500 40.75800 25.16600
N  -3.22300 34.33100 24.76300
C  -2.85500 35.30500 23.84100
C  -3.39300 34.95800 22.41400
C  -3.78700 33.44000 22.55100
C  -3.77600 33.28900 24.07200
C  -4.45800 35.91300 21.78300
C  -2.78800 32.46000 21.87700
C  -3.33300 31.25600 21.15100
N  -3.86100 32.67700 26.98500
C  -4.55900 31.95500 26.05900
C  -5.21600 30.76800 26.76200
C  -4.86600 30.86600 28.08700
C  -4.10200 32.07900 28.18800
C  -5.97600 29.69500 26.01600
C  -4.98700 30.30300 29.42000
O  -5.62100 29.26000 29.68400
C  -4.28200 31.35300 30.46900
C  -5.23900 32.03400 31.42500
O  -6.23700 32.64900 31.00000
O  -4.96700 31.73200 32.76000
C  -6.09800 32.02100 33.63400
H  -0.51900 37.13000 29.48800
H  -2.04400 37.20300 23.36800
H  -4.93400 31.57200 23.94000
H  -3.19000 33.89300 31.84400
H  -0.52800 34.94900 31.43900
H  -2.83400 36.13100 32.41600
H  -2.36000 36.92300 30.96100
H  -1.18900 36.74500 32.42000
H  -0.38600 32.98500 31.87000
H  -1.36700 31.74700 31.17900
H  -2.21200 31.11400 33.18900
H  -2.62300 32.76200 33.64700
H  -0.25500 40.16900 27.91900
H  1.01900 39.36100 27.27000
H  0.49900 38.75500 28.76900
H  0.66500 40.53600 25.48200
H  -0.91700 41.02300 26.06200
H  -0.46600 41.44500 24.32800
H  -2.50600 35.02500 21.78400
H  -4.81400 33.23900 22.24700
H  -4.50900 36.75000 22.47900
H  -5.42000 35.44200 21.58300
H  -4.08600 36.40600 20.88400
H  -2.24500 32.16900 22.77600
H  -2.06500 32.93000 21.21000
H  -2.75000 31.18100 20.23300
H  -4.37700 31.36900 20.85700
H  -3.03100 30.39100 21.74200
H  -5.97600 29.97600 24.96300
H  -7.01700 29.73400 26.33600
H  -5.54700 28.71900 26.24200
H  -3.48200 30.82800 30.99100
H  -6.59400 31.14900 34.05800
H  -6.86800 32.64200 33.17600
H  -5.69700 32.60500 34.46200


