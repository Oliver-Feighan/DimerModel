%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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


