%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
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
Mg 41.12800 41.85700 26.98600
C  40.12700 44.17900 29.64000
C  41.31000 39.41000 29.48700
C  42.36200 39.82100 24.75900
C  40.70600 44.36300 24.68300
N  40.96700 41.85000 29.29200
C  40.42100 42.87500 30.17300
C  40.33700 42.44400 31.65300
C  40.94600 40.93000 31.52100
C  41.11600 40.67200 30.01100
C  42.22400 40.71400 32.30800
C  38.87900 42.58800 32.18500
C  38.60800 41.92700 33.50300
H  37.80600 42.65600 34.55800
N  41.79000 39.86100 27.08100
C  41.78200 39.01000 28.19200
C  42.26600 37.68700 27.80200
C  42.58400 37.77600 26.38500
C  42.25600 39.19200 25.97700
C  42.31400 36.43000 28.74600
C  43.16500 36.78400 25.44300
O  43.35700 36.97800 24.25000
C  43.44600 35.41200 25.89800
N  41.39000 42.02100 25.04500
C  42.03300 41.09600 24.33000
C  42.18900 41.49200 22.84300
C  41.28100 42.70200 22.76900
C  41.15000 43.12300 24.28400
C  43.62900 41.80100 22.45000
C  39.92300 42.66800 22.11700
C  38.93300 41.74000 22.78800
N  40.64100 43.87900 27.09600
C  40.41100 44.75400 26.01900
C  40.01000 46.03800 26.53900
C  39.83600 45.76900 27.92100
C  40.24100 44.45600 28.22800
C  39.64600 47.28300 25.78400
C  39.34100 46.43200 29.16700
O  38.79000 47.54000 29.34100
C  39.55800 45.41800 30.28800
C  40.41800 46.03500 31.39300
O  41.61900 46.43600 31.36700
O  39.69400 45.95900 32.54600
C  40.38200 46.43100 33.76900
H  41.17000 38.61400 30.22000
H  42.83200 39.17800 24.01200
H  40.57600 45.21400 24.01100
H  40.93200 43.17900 32.19500
H  40.17200 40.21500 31.80000
H  42.06000 39.90200 33.01700
H  42.40500 41.63800 32.85600
H  43.16100 40.69100 31.75200
H  38.23700 42.28000 31.35900
H  38.65800 43.64900 32.30200
H  39.47300 41.65400 34.10700
H  37.99100 41.04400 33.33900
H  41.63300 35.71400 28.28500
H  41.97100 36.63500 29.76000
H  43.34800 36.08500 28.73700
H  43.92500 35.36100 26.87600
H  44.12500 35.02700 25.13800
H  42.45600 34.98600 26.05800
H  41.84700 40.64300 22.25100
H  41.73700 43.52300 22.21500
H  43.81400 41.38100 21.46200
H  44.33400 41.35700 23.15200
H  43.75500 42.87300 22.29400
H  40.06700 42.45700 21.05700
H  39.41200 43.62400 22.23000
H  38.39800 41.12300 22.06700
H  38.14200 42.27600 23.31300
H  39.38000 41.04100 23.49500
H  38.67000 47.63800 26.11400
H  39.60000 47.28100 24.69500
H  40.39300 48.03500 26.03900
H  38.59300 45.21900 30.75500
H  40.75500 45.57700 34.33400
H  39.56300 47.00000 34.20900
H  41.15300 47.17200 33.55400


