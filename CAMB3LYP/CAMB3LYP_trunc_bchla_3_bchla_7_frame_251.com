%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 1.27700 7.97600 26.81300
C  1.70500 10.25900 29.27500
C  1.82700 5.41400 29.09300
C  0.85500 5.69300 24.34900
C  1.21800 10.51800 24.46400
N  1.59300 7.92000 28.89400
C  1.77500 8.96200 29.74400
C  1.81600 8.52800 31.18600
C  2.13100 6.98700 31.08800
C  1.79700 6.74800 29.57900
C  3.53500 6.42400 31.56400
C  0.41400 8.75400 31.91200
C  0.55800 8.89600 33.45600
H  1.71800 8.24400 34.22700
N  1.40400 5.81900 26.67300
C  1.63700 4.96400 27.78200
C  1.59900 3.60000 27.34800
C  1.32800 3.64600 25.93400
C  1.17100 5.11200 25.62400
C  1.72200 2.42000 28.24700
C  1.18400 2.49900 24.94000
O  1.05900 2.62000 23.68600
C  1.20700 1.06700 25.44700
N  1.13600 8.07200 24.71300
C  0.77100 6.99400 23.90500
C  0.28500 7.38700 22.48000
C  0.38900 9.00800 22.61900
C  0.90800 9.24400 23.99000
C  1.23400 6.77600 21.44100
C  -0.88200 9.73900 22.30100
C  -0.72400 10.77600 21.22600
N  1.52400 9.91600 26.73700
C  1.45400 10.87900 25.73800
C  1.65900 12.16800 26.28700
C  1.81400 12.01500 27.69100
C  1.62000 10.64100 27.94600
C  1.67300 13.48200 25.60500
C  2.15100 12.69400 28.92300
O  2.45500 13.83200 29.22300
C  2.00700 11.52600 30.01000
C  3.06400 11.59000 31.05300
O  4.23800 11.79300 30.87600
O  2.59100 11.43200 32.29900
C  3.63300 11.47700 33.36400
H  2.18200 4.70500 29.84400
H  0.85400 4.95200 23.54700
H  1.09900 11.38200 23.80700
H  2.67100 8.94100 31.72100
H  1.35600 6.44600 31.63000
H  4.07900 7.25400 32.01500
H  4.16900 5.94900 30.81700
H  3.32200 5.63000 32.28000
H  -0.18400 7.89400 31.61000
H  -0.01700 9.66700 31.50100
H  -0.39200 8.49700 33.81200
H  0.49500 9.96300 33.66800
H  0.81200 2.12400 28.76900
H  2.52200 2.67200 28.94200
H  2.18600 1.56000 27.76500
H  2.26000 0.80200 25.54100
H  0.75700 0.48100 24.64600
H  0.71900 0.92900 26.41200
H  -0.69500 6.92500 22.36000
H  1.13100 9.36700 21.90500
H  2.13800 6.31200 21.83600
H  1.41900 7.51100 20.65800
H  0.69100 5.93800 21.00300
H  -1.27500 10.30600 23.14500
H  -1.71400 9.09900 22.00700
H  -1.48700 10.70900 20.45100
H  0.24000 10.68600 20.72600
H  -0.68400 11.71900 21.77200
H  2.02800 13.48300 24.57400
H  2.51200 13.99500 26.07600
H  0.69500 13.95600 25.68900
H  1.08300 11.88900 30.46100
H  3.95200 12.49700 33.58300
H  4.44500 10.77600 33.17000
H  3.18800 11.10700 34.28800
Mg 25.94600 0.51000 29.41700
C  27.75500 0.04800 32.42200
C  23.08400 0.98300 31.28700
C  24.18800 0.65700 26.52400
C  28.83500 -0.18800 27.55800
N  25.56000 0.39300 31.59400
C  26.46800 0.41400 32.65000
C  25.73500 0.61600 34.03900
C  24.24300 0.62700 33.61200
C  24.25900 0.61600 32.06700
C  23.30400 -0.50900 34.16600
C  26.25400 1.80200 34.95100
C  26.76400 1.39400 36.30000
H  26.45300 2.38700 37.48900
N  23.95800 0.76900 29.03000
C  22.92100 1.01800 29.90600
C  21.67800 1.14300 29.16100
C  21.94400 1.05800 27.76500
C  23.43600 0.83800 27.70800
C  20.32000 1.49500 29.80100
C  21.00100 1.19500 26.57400
O  21.36800 1.18900 25.41700
C  19.44600 1.46900 26.69600
N  26.43900 0.16900 27.31600
C  25.54800 0.29300 26.34100
C  26.18900 0.06900 24.92500
C  27.73100 -0.01000 25.26000
C  27.69500 0.01500 26.81800
C  25.54500 -1.13400 24.17700
C  28.51900 1.10000 24.59300
C  29.85400 0.64900 24.00300
N  27.89800 -0.11200 29.79400
C  28.99000 -0.26300 28.96100
C  30.08800 -0.69100 29.77200
C  29.69800 -0.55100 31.12700
C  28.32100 -0.19300 31.11800
C  31.42300 -1.10700 29.28200
C  30.17600 -0.59300 32.52000
O  31.23000 -0.96300 33.02900
C  28.81000 -0.12800 33.44600
C  28.50200 -1.15100 34.39600
O  27.57900 -1.93000 34.29600
O  29.31100 -1.04800 35.50500
C  29.06500 -2.12600 36.49500
H  22.18200 1.26000 31.83700
H  23.66500 0.64600 25.56600
H  29.70300 -0.30600 26.90500
H  25.82900 -0.35200 34.53000
H  23.76600 1.56200 33.90600
H  22.35300 -0.03400 34.40900
H  23.77900 -0.95200 35.04100
H  23.10600 -1.31400 33.45800
H  25.42400 2.50500 35.02300
H  27.09400 2.25600 34.42500
H  27.85000 1.29700 36.26800
H  26.30600 0.41800 36.45800
H  20.08400 2.49100 29.42800
H  20.46000 1.63600 30.87300
H  19.49000 0.83300 29.55100
H  19.31200 2.38400 27.27400
H  18.87900 0.57700 26.96400
H  19.10200 1.53100 25.66400
H  26.08400 0.92600 24.26000
H  28.04800 -1.00900 24.96000
H  26.34500 -1.78700 23.82700
H  24.86300 -0.83400 23.38100
H  24.93000 -1.70700 24.86900
H  28.70300 1.93600 25.26700
H  27.98200 1.61200 23.79300
H  30.10400 -0.36600 24.31200
H  30.60900 1.32100 24.41200
H  29.90400 0.79600 22.92400
H  31.58400 -1.10700 28.20400
H  31.68500 -2.00300 29.84600
H  32.04200 -0.24000 29.51400
H  29.03700 0.80200 33.96800
H  29.24700 -3.09800 36.03600
H  28.05200 -1.96800 36.86500
H  29.74500 -2.08500 37.34500


