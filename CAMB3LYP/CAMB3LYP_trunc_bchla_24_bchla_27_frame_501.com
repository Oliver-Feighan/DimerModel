%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
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
Mg -5.09800 24.63300 26.68400
C  -3.70200 26.36700 29.52200
C  -6.11300 22.14300 28.93500
C  -6.48000 23.09600 24.13300
C  -3.72900 27.06700 24.58100
N  -4.87100 24.33000 29.01500
C  -4.47600 25.23800 29.94600
C  -4.87000 24.81900 31.35300
C  -5.34300 23.26900 31.18300
C  -5.42800 23.22100 29.63000
C  -4.42000 22.22500 31.86000
C  -5.97200 25.74300 31.97300
C  -5.56400 26.39400 33.32000
H  -6.48200 25.91600 34.58700
N  -6.12400 22.92100 26.55900
C  -6.45400 21.95500 27.56900
C  -7.23800 20.80800 27.05300
C  -7.37200 21.15200 25.66800
C  -6.65700 22.41400 25.41400
C  -7.70400 19.61500 27.89000
C  -8.01000 20.37100 24.58300
O  -8.30600 20.80000 23.47300
C  -8.47200 18.93600 24.90400
N  -5.23600 25.12700 24.71500
C  -5.84300 24.30800 23.82600
C  -5.57400 24.69000 22.39700
C  -4.73400 25.98700 22.49800
C  -4.58300 26.13300 24.02900
C  -4.97900 23.58700 21.54700
C  -5.55400 27.25800 21.86300
C  -4.76800 28.10300 20.72000
N  -3.71900 26.23300 26.92800
C  -3.33300 27.18000 25.94900
C  -2.55400 28.22500 26.62300
C  -2.81300 27.99000 27.99900
C  -3.40400 26.73200 28.12600
C  -1.71900 29.35800 26.03400
C  -2.39200 28.46300 29.34800
O  -1.66000 29.39300 29.78700
C  -3.08100 27.46800 30.36000
C  -2.03500 26.94900 31.30200
O  -0.94600 26.57900 30.96500
O  -2.47500 27.07900 32.59000
C  -1.38400 26.91100 33.58300
H  -6.28600 21.33800 29.65200
H  -6.91400 22.63300 23.24400
H  -3.32400 27.85900 23.94800
H  -3.94900 24.82500 31.93600
H  -6.36100 23.10200 31.53600
H  -3.94100 21.58900 31.11500
H  -4.98800 21.61200 32.56100
H  -3.64700 22.62700 32.51500
H  -6.89500 25.19600 32.16700
H  -6.22100 26.57200 31.31000
H  -5.77000 27.46400 33.35800
H  -4.50500 26.32600 33.56700
H  -7.71300 19.74200 28.97300
H  -7.20000 18.67100 27.68300
H  -8.74200 19.45800 27.59500
H  -7.62000 18.31700 25.18500
H  -9.03200 18.48700 24.08300
H  -9.15600 18.98900 25.75100
H  -6.53500 24.86100 21.91000
H  -3.75100 25.85800 22.04600
H  -5.68000 23.16900 20.82400
H  -4.66700 22.69800 22.09500
H  -4.09700 23.98900 21.04800
H  -5.86900 27.99500 22.60100
H  -6.48500 26.95600 21.38400
H  -4.00100 27.42900 20.33800
H  -4.32800 29.01700 21.12100
H  -5.40200 28.28200 19.85200
H  -1.71900 30.24300 26.67000
H  -2.23000 29.70200 25.13400
H  -0.77900 28.87000 25.77600
H  -3.74400 28.08300 30.96800
H  -0.90300 25.95500 33.37700
H  -1.91800 26.79500 34.52600
H  -0.58100 27.64700 33.62200


