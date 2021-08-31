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
Mg -9.12600 18.73000 43.07000
C  -5.62300 18.44500 42.86700
C  -8.71100 22.08300 42.24600
C  -12.43100 19.03300 42.79300
C  -9.30300 15.26600 42.95700
N  -7.32000 20.02500 42.45300
C  -5.96700 19.76600 42.55400
C  -5.07000 20.98600 42.54900
C  -6.06600 22.15700 42.37300
C  -7.44800 21.42600 42.33500
C  -5.98600 23.25800 43.50700
C  -3.93600 21.03100 41.43000
C  -4.22600 20.35100 40.08200
H  -4.56600 21.25200 38.85800
N  -10.36200 20.30600 42.60500
C  -10.01500 21.60400 42.35900
C  -11.22400 22.42800 42.22000
C  -12.33600 21.54600 42.38200
C  -11.71100 20.20500 42.59800
C  -11.19700 23.93200 41.97500
C  -13.79400 22.10700 42.18800
O  -13.91800 23.28100 41.86800
C  -15.04600 21.14100 42.31600
N  -10.60200 17.36400 42.77700
C  -11.94900 17.71100 42.85200
C  -12.88700 16.37600 42.87100
C  -11.75200 15.22000 42.76300
C  -10.48500 16.00900 42.84500
C  -13.77000 16.28100 44.15000
C  -11.80700 14.33600 41.44400
C  -11.81900 12.80500 41.57400
N  -7.76200 17.11800 43.06200
C  -7.94400 15.76000 43.00400
C  -6.61700 15.11000 42.96500
C  -5.67900 16.17500 43.03100
C  -6.45600 17.36700 43.03800
C  -6.31500 13.64100 42.85700
C  -4.23100 16.45300 42.98800
O  -3.23100 15.72200 43.03800
C  -4.16400 17.97000 43.08300
C  -3.54200 18.42400 44.28100
O  -4.15400 18.45500 45.38000
O  -2.21800 18.84600 44.04800
C  -1.49600 19.49100 45.14000
H  -8.63900 23.16900 42.16100
H  -13.51700 19.11600 42.86600
H  -9.44500 14.18400 42.98300
H  -4.61500 20.94700 43.53800
H  -5.92300 22.50900 41.35200
H  -5.58200 24.18000 43.08900
H  -5.32600 22.85000 44.27300
H  -6.97200 23.48100 43.91600
H  -3.14100 20.46400 41.91600
H  -3.60700 22.06200 41.29800
H  -5.03700 19.62900 40.18300
H  -3.36900 19.72900 39.82700
H  -11.41500 24.23200 40.95000
H  -10.31200 24.38600 42.42100
H  -11.94400 24.43600 42.58800
H  -14.90900 20.56900 43.23400
H  -15.05000 20.49700 41.43700
H  -16.00600 21.63900 42.18000
H  -13.59700 16.50100 42.05400
H  -11.83300 14.61300 43.66500
H  -13.23300 15.66500 44.87100
H  -14.78200 15.93400 43.94100
H  -13.98700 17.24300 44.61400
H  -10.98000 14.66700 40.81500
H  -12.76100 14.64800 41.01900
H  -11.51100 12.31300 40.65100
H  -12.86600 12.60200 41.79800
H  -11.12600 12.60300 42.39100
H  -7.19000 13.22900 42.35500
H  -6.33400 13.32500 43.90000
H  -5.34100 13.39500 42.43300
H  -3.56700 18.21200 42.20400
H  -1.27600 18.88800 46.02100
H  -2.03900 20.37600 45.47300
H  -0.61300 19.94900 44.69400


