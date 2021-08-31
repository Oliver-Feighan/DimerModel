%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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
Mg 26.06500 50.95300 26.70600
C  23.95100 52.21500 29.08800
C  28.12100 49.83900 29.10000
C  28.15700 49.98600 24.24600
C  24.06600 52.63600 24.17800
N  26.14100 51.27800 28.89700
C  25.13500 51.71000 29.66700
C  25.49500 51.52800 31.23000
C  26.85400 50.74300 31.10400
C  27.11600 50.62500 29.58900
C  28.00100 51.35300 31.93300
C  24.40800 50.76600 32.04000
C  24.78700 50.21400 33.45300
H  24.26800 51.00700 34.71200
N  27.93000 50.07200 26.68200
C  28.56100 49.60000 27.73700
C  29.78100 48.94300 27.29300
C  29.90300 49.05100 25.87500
C  28.62800 49.72200 25.53900
C  30.86100 48.38000 28.20200
C  31.03700 48.64600 24.93900
O  31.09500 49.07700 23.79100
C  32.08400 47.71600 25.44400
N  26.12000 51.31300 24.51700
C  27.05500 50.76500 23.75000
C  26.75200 50.94300 22.25300
C  25.36900 51.58600 22.31600
C  25.20000 51.88900 23.75800
C  27.77300 51.70600 21.39600
C  24.29900 50.82000 21.58000
C  23.63800 51.55200 20.35800
N  24.37700 52.10200 26.60600
C  23.61600 52.65500 25.52200
C  22.39000 53.25700 25.99800
C  22.52400 53.06300 27.40100
C  23.71600 52.38900 27.69500
C  21.20200 53.78500 25.27400
C  21.72000 53.18300 28.59100
O  20.53300 53.46900 28.80800
C  22.59000 52.78500 29.75700
C  22.95900 54.06000 30.54700
O  23.29600 55.16300 30.09000
O  22.85900 53.76400 31.87400
C  23.13500 54.86100 32.85500
H  28.80700 49.32300 29.77500
H  28.85800 49.89200 23.41400
H  23.44600 52.99600 23.35500
H  25.65400 52.54800 31.57800
H  26.67500 49.72400 31.44700
H  28.52200 50.52300 32.40900
H  27.63000 51.99000 32.73600
H  28.68100 51.92100 31.29800
H  24.23600 49.84500 31.48300
H  23.53400 51.40200 32.18000
H  25.86900 50.31100 33.54800
H  24.40200 49.19900 33.54500
H  30.59600 48.49500 29.25300
H  31.79900 48.90500 28.02300
H  31.04200 47.32300 28.01100
H  31.53500 46.89000 25.89700
H  32.83600 48.00900 26.17600
H  32.72300 47.32900 24.65100
H  26.66800 49.94900 21.81500
H  25.43600 52.54200 21.79700
H  27.79000 51.15500 20.45600
H  28.73100 51.77600 21.91100
H  27.54500 52.76300 21.25900
H  23.50000 50.72800 22.31500
H  24.63300 49.85300 21.20500
H  24.06100 52.55100 20.24900
H  22.57600 51.65700 20.58200
H  23.65600 51.00200 19.41700
H  21.27300 53.54700 24.21300
H  21.29800 54.86900 25.32400
H  20.34800 53.31100 25.76000
H  21.96100 52.11200 30.34000
H  24.03500 55.44400 32.65800
H  23.15200 54.54000 33.89600
H  22.32000 55.57300 32.72600


