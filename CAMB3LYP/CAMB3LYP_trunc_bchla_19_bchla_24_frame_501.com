%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
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
Mg 0.26000 43.34200 24.75600
C  2.25300 43.09900 27.61800
C  -2.44200 42.42000 26.69400
C  -1.57400 42.87700 21.91500
C  3.14600 43.60500 22.86000
N  -0.01100 42.92700 26.88700
C  0.88400 43.10000 27.86800
C  0.29200 43.09400 29.28500
C  -1.20300 42.68500 28.99500
C  -1.29800 42.75900 27.37000
C  -1.68200 41.32300 29.54200
C  0.44600 44.41100 30.04600
C  0.46000 44.23200 31.57500
H  0.30000 42.82500 32.16800
N  -1.74500 42.84700 24.34000
C  -2.74500 42.55800 25.29700
C  -3.99500 42.34000 24.59700
C  -3.74800 42.40500 23.17500
C  -2.33100 42.75300 23.08200
C  -5.28700 42.18000 25.35000
C  -4.75800 42.32200 21.94400
O  -4.44000 42.53400 20.79900
C  -6.18400 41.93600 22.27200
N  0.72800 43.16000 22.71000
C  -0.18900 43.01100 21.71100
C  0.47200 43.04900 20.28400
C  1.92800 43.34600 20.63900
C  1.92700 43.36400 22.15600
C  0.26000 41.87700 19.24300
C  2.38800 44.69400 20.03600
C  2.04300 46.00200 20.78200
N  2.32100 43.50900 25.09100
C  3.38200 43.61100 24.22900
C  4.61600 43.69300 24.93600
C  4.22700 43.48400 26.33500
C  2.89200 43.36300 26.33900
C  5.99200 43.88800 24.47400
C  4.71200 43.25500 27.68900
O  5.85600 43.34100 28.14800
C  3.47500 43.00100 28.60000
C  3.56600 41.71300 29.28000
O  3.33400 40.61600 28.78300
O  4.17100 41.89000 30.49900
C  4.75700 40.73100 31.20000
H  -3.23100 42.10200 27.37900
H  -2.18100 42.91300 21.00800
H  4.04900 43.78500 22.27300
H  0.67700 42.26400 29.87700
H  -1.92300 43.45600 29.26700
H  -2.57300 41.63600 30.08600
H  -0.94100 40.89700 30.21800
H  -1.87200 40.62800 28.72400
H  -0.32500 45.13200 29.77700
H  1.38300 44.89600 29.77100
H  -0.37400 44.87500 31.85700
H  1.32100 44.72000 32.03400
H  -5.00100 42.06900 26.39600
H  -5.81000 41.30000 24.97500
H  -5.96200 43.02800 25.22500
H  -6.37400 40.88800 22.50600
H  -6.81500 42.27400 21.45000
H  -6.38700 42.55900 23.14300
H  -0.03500 43.91100 19.85100
H  2.51500 42.52100 20.23700
H  1.17400 41.51200 18.77600
H  -0.34600 42.26100 18.42300
H  -0.30900 41.09200 19.74200
H  1.97300 44.71300 19.02900
H  3.46900 44.59400 19.93400
H  1.75200 46.69500 19.99300
H  2.89700 46.44600 21.29300
H  1.29400 45.88900 21.56500
H  6.47500 44.58100 25.16200
H  6.09900 44.26000 23.45500
H  6.57000 42.97800 24.63600
H  3.48200 43.84300 29.29200
H  4.58400 39.79600 30.66600
H  4.26600 40.64900 32.17000
H  5.80700 41.02000 31.22900


