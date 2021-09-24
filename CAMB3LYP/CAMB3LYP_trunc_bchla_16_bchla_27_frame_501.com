%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 40.80400 41.48400 27.15800
C  39.93800 43.96600 29.51800
C  41.67100 39.43500 29.59200
C  42.41200 39.59100 24.81200
C  40.80000 44.20800 24.65700
N  40.84000 41.74700 29.31900
C  40.23400 42.70800 30.04600
C  40.11400 42.23500 31.50800
C  40.62500 40.83000 31.52200
C  41.12800 40.63200 30.06000
C  41.62900 40.55000 32.71400
C  38.74000 42.40700 32.23800
C  38.20000 41.35400 33.29800
H  37.58800 41.87900 34.60500
N  41.67300 39.62800 27.19200
C  41.94100 38.92200 28.34300
C  42.50600 37.66000 27.95800
C  42.79100 37.65600 26.53700
C  42.34000 39.00600 26.13300
C  42.84100 36.69600 29.02300
C  43.34600 36.53900 25.61900
O  43.54900 36.75600 24.47000
C  43.39500 35.12400 26.09600
N  41.54800 41.90500 25.04000
C  42.04300 40.86300 24.33500
C  42.34300 41.24800 22.85100
C  41.79900 42.70700 22.76800
C  41.36100 43.00200 24.22400
C  43.78600 41.15300 22.43200
C  40.63300 43.00600 21.76200
C  39.24100 42.61000 22.12200
N  40.45200 43.63600 26.97900
C  40.44000 44.51600 25.93500
C  39.94600 45.79200 26.39200
C  39.72700 45.66500 27.78900
C  40.09700 44.34100 28.10100
C  39.86300 47.06500 25.62600
C  39.33400 46.33800 28.99400
O  38.92600 47.46900 29.23200
C  39.40400 45.25700 30.16400
C  40.16300 45.75100 31.38400
O  41.33900 46.00900 31.38400
O  39.30000 45.80400 32.46700
C  39.93400 46.30500 33.68700
H  41.79100 38.62500 30.31500
H  42.89000 38.90600 24.10900
H  40.69600 44.98800 23.90000
H  40.78700 42.96300 31.96100
H  39.78000 40.17800 31.74500
H  41.64600 41.41100 33.38200
H  42.57700 40.23500 32.27800
H  41.25100 39.76900 33.37400
H  37.98000 42.45900 31.45800
H  38.66600 43.37100 32.74200
H  39.06500 40.72800 33.51800
H  37.46000 40.71500 32.81600
H  43.75800 36.20700 28.69500
H  42.05100 35.94600 29.03300
H  42.98100 37.00500 30.05900
H  44.39900 35.05700 26.51300
H  43.13900 34.46100 25.26900
H  42.74200 34.89700 26.93900
H  41.71500 40.69400 22.15300
H  42.49200 43.53200 22.60600
H  44.38600 40.63600 23.18200
H  44.24800 42.13400 22.31800
H  43.78100 40.59500 21.49600
H  40.84900 42.53700 20.80200
H  40.70100 44.05800 21.48400
H  39.26600 41.80800 22.86000
H  38.56300 42.19300 21.37700
H  38.76600 43.53300 22.45600
H  40.81300 47.42600 25.23300
H  39.70300 47.84000 26.37600
H  39.09100 47.02600 24.85700
H  38.36800 45.04000 30.42500
H  39.84500 45.31900 34.14200
H  39.21300 47.03100 34.06200
H  40.94500 46.67900 33.52400
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


