%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg -1.96800 16.99100 27.07800
C  -2.18500 14.77300 29.87900
C  -3.22300 19.43300 29.10400
C  -2.21300 18.84100 24.51400
C  -1.82200 14.16700 25.04100
N  -2.58500 16.98200 29.20400
C  -2.56200 16.08300 30.18300
C  -2.90400 16.68200 31.58000
C  -3.74100 17.95500 31.11000
C  -3.12300 18.15900 29.73800
C  -5.23100 17.67600 30.94100
C  -1.62300 17.10100 32.37700
C  -1.84400 17.38800 33.86500
H  -0.80000 17.02400 34.81100
N  -2.15300 18.95200 26.92600
C  -2.67100 19.79100 27.85500
C  -2.58800 21.14700 27.27100
C  -2.15500 21.00600 25.85300
C  -2.18200 19.57400 25.69800
C  -2.85500 22.39700 28.05200
C  -1.78700 22.00400 24.80100
O  -1.66500 21.68100 23.57000
C  -1.65800 23.40500 25.24100
N  -2.23300 16.57800 25.04900
C  -2.24700 17.56000 24.14100
C  -2.16400 17.01500 22.74600
C  -1.84400 15.45800 22.89800
C  -1.88200 15.35900 24.44300
C  -3.41900 17.19600 21.88100
C  -0.59500 14.92100 22.15500
C  0.65800 15.50100 22.77200
N  -1.81800 14.84700 27.38100
C  -1.79500 13.85500 26.39000
C  -1.79500 12.55300 27.08200
C  -1.83500 12.88000 28.43600
C  -1.88500 14.29600 28.58000
C  -1.81900 11.27600 26.45500
C  -1.87200 12.36200 29.79400
O  -1.76000 11.24400 30.29900
C  -2.22300 13.56300 30.77900
C  -1.25100 13.50900 31.92600
O  -0.03700 13.73600 31.74700
O  -1.93000 13.43500 33.08100
C  -1.02500 13.40800 34.26100
H  -3.82700 20.21100 29.57500
H  -2.18600 19.44400 23.60400
H  -1.69200 13.24400 24.47300
H  -3.52200 16.04200 32.21000
H  -3.55300 18.82200 31.74300
H  -5.45200 16.61400 31.04700
H  -5.47900 18.01800 29.93600
H  -5.82000 18.29200 31.62100
H  -1.29900 18.08700 32.04400
H  -0.83800 16.35200 32.26700
H  -2.81200 17.01700 34.20300
H  -1.90600 18.46800 33.99900
H  -3.31200 22.17700 29.01600
H  -3.61900 22.95400 27.51000
H  -1.96100 23.01100 28.16200
H  -2.60300 23.73500 25.67300
H  -1.30300 23.96100 24.37400
H  -0.92000 23.45700 26.04200
H  -1.36000 17.49200 22.18700
H  -2.66000 14.83500 22.53100
H  -3.23500 17.73100 20.94900
H  -4.17200 17.70500 22.48200
H  -3.93500 16.27100 21.62500
H  -0.61200 15.23200 21.11100
H  -0.51600 13.83400 22.17600
H  0.74600 16.46900 22.28000
H  1.51900 14.87000 22.54900
H  0.52700 15.68800 23.83800
H  -1.10100 10.65800 26.99400
H  -1.56400 11.34900 25.39800
H  -2.85900 10.99200 26.61500
H  -3.25100 13.41000 31.10600
H  -1.27000 12.51300 34.83300
H  -1.13500 14.37200 34.75800
H  0.03200 13.23700 34.05400
Mg 17.05400 -2.28800 27.75000
C  16.54800 -0.35400 30.72600
C  18.85200 -4.52900 29.65700
C  18.01300 -3.67000 24.94900
C  15.39200 0.22200 26.03300
N  17.88600 -2.27700 29.96200
C  17.41200 -1.43700 30.93300
C  17.78200 -2.00100 32.33100
C  18.58500 -3.32700 31.99600
C  18.48400 -3.42400 30.43000
C  20.10200 -3.41100 32.45300
C  16.57800 -2.28600 33.26800
C  16.68400 -1.85000 34.68600
H  18.05300 -1.52300 35.39600
N  18.30500 -3.91800 27.31000
C  18.85900 -4.81900 28.21900
C  19.64100 -5.82500 27.47500
C  19.29300 -5.63200 26.08300
C  18.51600 -4.37100 26.05800
C  20.53100 -6.87100 28.07500
C  19.63900 -6.46300 24.82500
O  19.04500 -6.38400 23.72900
C  20.63000 -7.54600 24.97700
N  16.79100 -1.75400 25.84400
C  17.37300 -2.49900 24.84200
C  16.84500 -2.07000 23.44500
C  15.70300 -1.04100 23.86200
C  15.95300 -0.84300 25.34700
C  17.88000 -1.47900 22.40300
C  14.26500 -1.46100 23.60200
C  13.50200 -0.61900 22.55800
N  15.98100 -0.53800 28.25500
C  15.38600 0.39400 27.43800
C  14.87800 1.51400 28.24100
C  15.27300 1.20700 29.53100
C  15.93800 -0.03900 29.49700
C  14.17200 2.75100 27.79300
C  15.18900 1.69900 30.85900
O  14.57400 2.68200 31.36200
C  15.98100 0.75800 31.66400
C  17.01200 1.50300 32.46800
O  17.97500 2.03100 32.02700
O  16.72000 1.42300 33.76600
C  17.36800 2.29700 34.68400
H  19.44500 -5.28100 30.18200
H  18.21500 -4.04600 23.94300
H  14.94000 1.04300 25.47400
H  18.42200 -1.24900 32.79400
H  18.07000 -4.18900 32.41900
H  20.43900 -2.60400 33.10300
H  20.70600 -3.31100 31.55000
H  20.31100 -4.38200 32.90300
H  16.34900 -3.34900 33.35200
H  15.66000 -1.80400 32.93300
H  16.30900 -2.63700 35.34100
H  16.04900 -0.97200 34.80100
H  21.52400 -6.79800 27.63200
H  20.09600 -7.83300 27.80400
H  20.47400 -6.87000 29.16400
H  21.57900 -7.08500 25.25100
H  20.75500 -7.96300 23.97800
H  20.26200 -8.35300 25.61000
H  16.37000 -2.96000 23.03100
H  15.88800 -0.12500 23.30100
H  18.88000 -1.46500 22.83500
H  17.57100 -0.51800 21.99000
H  17.86400 -2.13300 21.53100
H  13.63500 -1.61600 24.47800
H  14.30700 -2.45900 23.16700
H  12.47900 -0.33500 22.80400
H  13.29800 -1.28100 21.71600
H  13.99200 0.30400 22.25100
H  13.75500 2.65600 26.79000
H  14.88900 3.56000 27.65400
H  13.49700 3.03000 28.60200
H  15.27200 0.27100 32.33300
H  16.91900 3.28000 34.53500
H  18.41700 2.43200 34.42000
H  17.17400 2.00000 35.71500


