%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 1.29600 7.71300 26.88100
C  1.46900 9.82700 29.67000
C  2.42600 5.06900 28.94700
C  1.47300 5.72500 24.13400
C  0.92500 10.50000 24.87700
N  1.69900 7.49100 29.04900
C  1.51500 8.49600 30.00300
C  1.43800 7.84300 31.43000
C  1.83300 6.38800 31.15200
C  2.05200 6.30500 29.60100
C  3.12000 6.01700 31.93700
C  0.04000 8.06800 31.99300
C  -0.12000 8.13500 33.50300
H  1.05800 7.66300 34.38500
N  1.84800 5.68100 26.53500
C  2.31100 4.77600 27.50900
C  2.63500 3.49600 26.90400
C  2.31300 3.59100 25.49300
C  1.87500 5.04400 25.33500
C  3.03900 2.24400 27.69300
C  2.37500 2.55400 24.38800
O  2.03400 2.84400 23.24900
C  2.64700 1.14500 24.75200
N  1.13700 8.09300 24.78900
C  1.15400 7.08500 23.87300
C  0.77800 7.64200 22.47500
C  0.53200 9.18300 22.70900
C  0.76500 9.29400 24.23800
C  1.89300 7.33600 21.45900
C  -0.89600 9.67800 22.35700
C  -0.86300 10.38200 21.00400
N  1.30900 9.79300 27.10700
C  1.22300 10.82400 26.23200
C  1.50800 12.06900 26.98600
C  1.50300 11.75700 28.31400
C  1.42300 10.34100 28.32600
C  1.55300 13.42200 26.29600
C  1.55800 12.26300 29.69200
O  1.68500 13.41500 30.08900
C  1.44800 11.03500 30.57500
C  2.56000 11.13900 31.53300
O  3.76600 10.89700 31.38900
O  2.04800 11.65700 32.70500
C  2.99200 11.90900 33.77200
H  2.69600 4.26200 29.63100
H  1.45200 4.99400 23.32300
H  0.89300 11.33000 24.16800
H  2.21800 8.40600 31.94400
H  1.06000 5.77700 31.61800
H  4.00400 5.85500 31.32000
H  2.97100 4.99800 32.29300
H  3.38200 6.77100 32.67900
H  -0.56100 7.28300 31.53300
H  -0.37400 8.97400 31.55200
H  -1.02500 7.59600 33.78000
H  -0.31800 9.19800 33.64500
H  3.22100 2.42800 28.75200
H  3.95300 1.82100 27.27400
H  2.28900 1.47200 27.52500
H  1.81900 0.77800 25.35900
H  3.66100 1.02700 25.13400
H  2.61600 0.48900 23.88200
H  -0.10600 7.10300 22.13500
H  1.34500 9.68200 22.18000
H  2.33200 8.19000 20.94400
H  1.31000 6.77100 20.73200
H  2.77000 6.84600 21.88100
H  -1.17200 10.51700 22.99600
H  -1.64500 8.88600 22.34300
H  -1.71200 9.95300 20.47200
H  -0.02300 10.09700 20.37100
H  -0.94100 11.46900 21.03900
H  1.59500 13.38900 25.20700
H  2.38500 14.03800 26.63500
H  0.58900 13.91000 26.44400
H  0.47700 11.06900 31.06800
H  3.81000 12.49700 33.35500
H  3.40600 11.00100 34.21100
H  2.51900 12.68100 34.37800
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


