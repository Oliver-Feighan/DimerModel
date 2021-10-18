%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 8.37300 3.02000 27.58300
C  9.66800 1.59900 30.62700
C  6.70500 5.26800 29.53700
C  6.93600 4.14100 24.76700
C  9.73100 0.40000 25.89700
N  8.26000 3.41400 29.92600
C  8.92300 2.70600 30.88600
C  8.72800 3.42100 32.16900
C  7.67000 4.51900 31.87000
C  7.43900 4.42400 30.33700
C  6.33800 4.27100 32.64700
C  9.97700 3.96400 32.81700
C  9.90300 4.49900 34.23800
H  10.91500 4.00800 35.18600
N  7.11500 4.54900 27.20100
C  6.43400 5.32900 28.17600
C  5.69100 6.36700 27.43700
C  5.70300 5.96800 26.10300
C  6.57300 4.78400 25.98300
C  4.94900 7.46900 28.16400
C  5.04400 6.62700 24.87200
O  4.92600 6.07300 23.76300
C  4.55800 8.00900 25.03700
N  8.26900 2.33700 25.58000
C  7.65300 2.93000 24.57800
C  7.79400 2.17100 23.23000
C  8.87800 1.15300 23.59100
C  9.00200 1.33000 25.08700
C  6.44500 1.50000 22.80300
C  10.20800 1.29100 22.88500
C  10.86600 -0.01000 22.39200
N  9.53800 1.41400 28.13000
C  10.02500 0.42800 27.30000
C  10.88100 -0.54000 28.05000
C  10.69600 -0.06300 29.33200
C  9.85900 1.05900 29.36100
C  11.78700 -1.59900 27.50900
C  11.07900 -0.33300 30.71000
O  11.83900 -1.20400 31.11800
C  10.43500 0.73200 31.62700
C  9.77600 0.01400 32.68800
O  8.68700 -0.60500 32.61700
O  10.41500 0.23000 33.87100
C  9.76200 -0.17400 35.09700
H  6.23400 6.04600 30.14100
H  6.34700 4.39700 23.88400
H  9.91800 -0.50400 25.31400
H  8.28800 2.65400 32.80500
H  8.01700 5.51600 32.14200
H  5.84500 5.03100 33.25400
H  6.57700 3.45100 33.32400
H  5.53900 3.81600 32.06100
H  10.35300 4.83500 32.28000
H  10.84000 3.32300 32.63900
H  8.94300 4.25600 34.69100
H  10.21000 5.54200 34.16500
H  3.93400 7.63800 27.80500
H  5.45900 8.43200 28.14900
H  4.76300 7.27700 29.22100
H  5.36300 8.60800 25.46300
H  3.75000 7.88400 25.75700
H  4.37900 8.39200 24.03300
H  8.13900 2.90800 22.50400
H  8.61200 0.10400 23.46400
H  6.50600 0.41300 22.84700
H  6.19800 1.73200 21.76700
H  5.70700 1.85800 23.52100
H  10.93600 1.82600 23.49400
H  10.03500 2.08200 22.15600
H  11.47800 0.25700 21.53100
H  10.15100 -0.81800 22.23400
H  11.45900 -0.28500 23.26400
H  12.83900 -1.31500 27.47600
H  11.33000 -1.91300 26.57100
H  11.63300 -2.52100 28.07100
H  11.28200 1.20800 32.12100
H  8.73800 -0.52700 34.97200
H  9.83700 0.61100 35.84900
H  10.27900 -1.00300 35.58100
Mg 16.95800 -2.04300 28.05800
C  15.67900 0.03700 30.54100
C  18.78700 -3.72100 30.44000
C  18.31300 -3.78900 25.50500
C  15.32700 0.02400 25.64200
N  17.25900 -1.76300 30.19400
C  16.62800 -0.88400 31.01900
C  16.98400 -1.14100 32.48900
C  18.10000 -2.26900 32.35300
C  18.02400 -2.61400 30.86200
C  19.48200 -1.75800 32.74200
C  15.83300 -1.72700 33.38200
C  15.87500 -1.41800 34.85600
H  17.10800 -0.90300 35.45200
N  18.26300 -3.57700 28.00200
C  18.92400 -4.12800 29.11800
C  19.62000 -5.29300 28.63800
C  19.47300 -5.38800 27.23100
C  18.65200 -4.21500 26.84200
C  20.55500 -6.11300 29.50100
C  19.94500 -6.39800 26.27100
O  19.91000 -6.24200 25.06900
C  20.68300 -7.69100 26.77500
N  16.72600 -2.04100 25.82800
C  17.40600 -2.83200 25.05300
C  17.12300 -2.53800 23.54900
C  16.12900 -1.39100 23.56700
C  16.05500 -1.03300 25.09500
C  18.49500 -2.31000 22.75600
C  14.77700 -1.75000 22.91500
C  14.57500 -1.34500 21.40700
N  15.59300 -0.54700 28.03400
C  15.01200 0.20800 26.98500
C  14.13500 1.19500 27.48900
C  14.45500 1.19400 28.91300
C  15.31700 0.08300 29.13900
C  13.41400 2.25000 26.67000
C  14.16900 1.92500 30.09200
O  13.43900 2.88400 30.30900
C  14.95100 1.07800 31.29700
C  15.78000 2.10000 32.01800
O  16.65500 2.79200 31.48700
O  15.36200 2.29400 33.36800
C  16.13600 3.30900 34.14900
H  19.34000 -4.20600 31.24800
H  18.93000 -4.26900 24.74300
H  14.76600 0.65500 24.95000
H  17.39400 -0.20500 32.86900
H  17.90900 -3.13400 32.98800
H  19.41600 -0.81500 33.28400
H  20.16500 -1.70100 31.89400
H  19.97000 -2.35300 33.51400
H  15.82100 -2.80600 33.23200
H  14.89600 -1.40800 32.92500
H  15.72700 -2.36200 35.38000
H  15.08500 -0.70500 35.09200
H  21.59900 -6.03400 29.19800
H  20.20800 -7.14700 29.48400
H  20.60300 -5.83600 30.55400
H  21.71600 -7.39100 26.95300
H  20.76400 -8.41600 25.96500
H  20.25400 -8.04400 27.71200
H  16.67800 -3.47000 23.20000
H  16.55100 -0.48300 23.13900
H  18.49400 -1.34000 22.25900
H  18.63800 -3.00100 21.92500
H  19.26300 -2.34000 23.53000
H  14.01500 -1.24300 23.50800
H  14.64200 -2.82300 23.05200
H  14.22000 -0.31400 21.38800
H  13.95900 -1.97800 20.76800
H  15.61200 -1.30800 21.07500
H  12.51000 1.73500 26.34600
H  13.99600 2.54300 25.79700
H  13.31000 3.23400 27.12900
H  14.20300 0.59900 31.93000
H  16.03100 4.30700 33.72500
H  17.20500 3.14900 34.29200
H  15.57400 3.26700 35.08200


