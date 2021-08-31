%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
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
Mg 9.20100 48.63300 24.80700
C  7.00600 48.65300 27.54500
C  11.64200 49.94100 26.66500
C  11.05000 48.61100 22.15100
C  6.24200 48.21900 22.76700
N  9.22500 49.31100 26.85200
C  8.33600 49.09200 27.82900
C  8.78900 49.67200 29.14400
C  10.10600 50.44700 28.73700
C  10.36100 49.87000 27.35800
C  10.08400 52.01000 28.74700
C  9.14700 48.49500 30.18800
C  8.51500 48.70500 31.60300
H  8.63500 47.36700 32.51500
N  11.09800 48.95000 24.52400
C  12.00800 49.46400 25.43200
C  13.32900 49.31600 24.88600
C  13.24200 48.87600 23.50200
C  11.75200 48.67100 23.33700
C  14.54800 49.78900 25.72900
C  14.23500 48.59100 22.38100
O  13.96800 48.12400 21.25900
C  15.64500 48.82300 22.79500
N  8.62300 48.59400 22.78300
C  9.67100 48.60700 21.89600
C  9.16600 48.63700 20.45300
C  7.69000 48.19700 20.63800
C  7.46000 48.39900 22.13700
C  9.46200 50.04800 19.80600
C  7.33800 46.72100 20.17100
C  5.87200 46.44300 19.80200
N  7.05000 48.42700 25.05400
C  6.04600 48.27100 24.14700
C  4.77900 48.04000 24.80400
C  5.15700 48.11700 26.15700
C  6.52300 48.41500 26.26400
C  3.35900 47.77300 24.25400
C  4.63600 47.99600 27.49000
O  3.59000 47.60400 27.93700
C  5.82500 48.51800 28.45200
C  5.77300 47.51400 29.51800
O  6.28400 46.36600 29.55000
O  5.03000 48.07600 30.48700
C  4.75700 47.29300 31.69300
H  12.47200 50.28300 27.28700
H  11.67900 48.52800 21.26200
H  5.39900 48.00100 22.10800
H  8.05900 50.30600 29.64800
H  11.00500 50.19100 29.29800
H  10.13500 52.34800 27.71200
H  10.98900 52.30100 29.28000
H  9.22300 52.51800 29.18000
H  10.20800 48.27500 30.31200
H  8.82700 47.56600 29.71600
H  7.46000 48.96400 31.51600
H  9.00400 49.45600 32.22300
H  14.89000 48.94400 26.32500
H  14.28100 50.45800 26.54800
H  15.30300 50.38600 25.21800
H  15.98500 48.19800 23.62000
H  15.64000 49.89000 23.02100
H  16.36500 48.70300 21.98500
H  9.75900 47.88800 19.92800
H  7.04700 48.83900 20.03600
H  8.53800 50.61700 19.70400
H  9.90700 50.08500 18.81200
H  10.06100 50.74200 20.39600
H  7.49400 46.14100 21.08100
H  7.92900 46.29300 19.36100
H  5.39900 45.74600 20.49500
H  5.63800 46.05100 18.81200
H  5.34900 47.38100 19.98300
H  3.34300 48.21100 23.25600
H  2.59900 48.22600 24.89200
H  3.25800 46.69000 24.18400
H  5.55100 49.48200 28.88100
H  5.30400 47.73800 32.52400
H  4.92100 46.22700 31.53900
H  3.71000 47.43600 31.96000


