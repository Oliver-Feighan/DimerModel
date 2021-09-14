%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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
Mg 25.47400 50.80800 26.36000
C  23.42800 51.77500 29.14000
C  27.78000 49.80000 28.46900
C  27.48600 50.74800 23.72700
C  22.89100 52.35600 24.24800
N  25.55500 50.75400 28.55700
C  24.64000 51.22200 29.52100
C  25.12200 51.01800 31.00700
C  26.56400 50.28500 30.68300
C  26.67300 50.24400 29.16200
C  27.79400 51.06900 31.31300
C  24.14500 50.22000 31.93900
C  24.75600 49.66700 33.22100
H  24.21600 50.30100 34.51600
N  27.35300 50.24000 26.12300
C  28.22400 49.80700 27.09500
C  29.50900 49.42000 26.55200
C  29.38500 49.70200 25.13400
C  27.97500 50.19600 24.94200
C  30.61900 48.81900 27.27800
C  30.38700 49.50900 24.00900
O  30.13300 49.66800 22.84800
C  31.78500 49.18000 24.46100
N  25.24600 51.52100 24.28400
C  26.25500 51.32600 23.39800
C  25.83300 51.76400 21.97000
C  24.32300 52.13900 22.16700
C  24.10500 51.95100 23.65600
C  26.64800 52.90800 21.27600
C  23.28700 51.29000 21.34900
C  23.20200 51.69600 19.86100
N  23.59000 51.86200 26.56100
C  22.66400 52.33500 25.66000
C  21.52300 52.77100 26.36200
C  21.75600 52.64400 27.68000
C  23.06700 52.09200 27.80400
C  20.26700 53.21900 25.62600
C  21.30400 52.86100 29.05800
O  20.28500 53.34300 29.48100
C  22.42600 52.30600 30.05900
C  23.02400 53.48800 30.75300
O  23.93000 54.23000 30.44700
O  22.33300 53.57900 31.94500
C  22.72800 54.66700 32.87100
H  28.62400 49.48900 29.08800
H  28.23800 50.84300 22.94000
H  22.18300 52.67600 23.48000
H  25.36700 52.00000 31.41200
H  26.54900 49.29900 31.14800
H  28.47500 50.35700 31.77900
H  27.60200 51.88000 32.01500
H  28.41300 51.49500 30.52400
H  23.78800 49.47000 31.23200
H  23.21100 50.69500 32.24000
H  25.83400 49.82400 33.19700
H  24.45000 48.62700 33.33100
H  30.63000 47.83900 26.80100
H  30.42700 48.80700 28.35100
H  31.61000 49.25400 27.15000
H  31.91600 48.09800 24.46700
H  32.15000 49.82600 25.26000
H  32.40300 49.70200 23.73100
H  25.75100 50.83500 21.40500
H  24.20500 53.21400 22.02900
H  26.08400 53.79400 20.98300
H  26.95200 52.49600 20.31400
H  27.50400 53.28000 21.83900
H  22.29300 51.43100 21.77200
H  23.59200 50.25400 21.50300
H  23.57200 52.71200 19.72200
H  22.15600 51.57400 19.58000
H  23.82700 51.09700 19.19800
H  19.37700 53.06800 26.23700
H  20.20500 52.68100 24.68000
H  20.36300 54.25800 25.31000
H  22.01700 51.57800 30.76000
H  21.79500 55.04700 33.28700
H  23.27800 55.48100 32.39900
H  23.31400 54.08500 33.58300


