%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 26.07000 0.32300 29.17500
C  27.93000 0.00300 32.17200
C  23.23700 0.37300 31.22200
C  24.20000 0.73500 26.50800
C  28.89900 -0.55900 27.35600
N  25.68800 0.23300 31.42600
C  26.58700 0.25300 32.44900
C  25.99300 0.25100 33.86800
C  24.43500 0.24800 33.45400
C  24.48700 0.27000 31.95500
C  23.63400 -1.02300 33.94300
C  26.39900 1.50900 34.70300
C  27.30600 1.47800 35.92700
H  27.07500 2.41400 37.10600
N  24.00400 0.60600 28.87400
C  22.98600 0.52400 29.83500
C  21.72400 0.61300 29.13500
C  21.97300 0.88000 27.79100
C  23.44300 0.80400 27.68300
C  20.46800 0.64800 29.93700
C  21.02800 1.08200 26.60700
O  21.39200 1.28100 25.44400
C  19.58500 1.18400 26.91000
N  26.42900 0.02500 27.20000
C  25.57000 0.46200 26.23000
C  26.30100 0.60400 24.85200
C  27.70600 0.10100 25.15700
C  27.68000 -0.24400 26.67200
C  25.55900 -0.10800 23.67400
C  28.73300 1.18300 24.79000
C  29.96300 0.77300 24.06500
N  28.00500 -0.28200 29.61000
C  29.05600 -0.58000 28.78000
C  30.24500 -0.86800 29.52500
C  29.92400 -0.56000 30.93700
C  28.52000 -0.19800 30.89800
C  31.58600 -1.33700 28.93500
C  30.40700 -0.48700 32.33900
O  31.55800 -0.58600 32.79700
C  29.11500 -0.11100 33.17500
C  29.05800 -1.28400 34.10300
O  28.53100 -2.37500 33.92200
O  29.56000 -0.97800 35.26800
C  29.75400 -2.03800 36.28700
H  22.46700 0.52900 31.98000
H  23.69600 0.81600 25.54300
H  29.79300 -0.76400 26.76400
H  26.26300 -0.69600 34.33700
H  23.82700 1.11200 33.72400
H  24.21800 -1.88400 34.26700
H  22.99700 -1.39000 33.13900
H  22.99100 -0.72400 34.77100
H  25.46900 1.94600 35.06700
H  26.90100 2.23700 34.06700
H  28.23400 1.84700 35.49000
H  27.27900 0.42100 36.19200
H  19.92700 1.58400 29.80500
H  20.67800 0.53600 31.00100
H  19.76300 -0.13800 29.66800
H  19.43300 1.90200 27.71600
H  19.18900 0.27500 27.36300
H  18.97400 1.53500 26.07900
H  26.25300 1.68100 24.69200
H  27.96600 -0.74800 24.52500
H  24.56100 -0.46700 23.92500
H  26.08600 -0.97800 23.28400
H  25.56500 0.53500 22.79400
H  29.00400 1.58600 25.76500
H  28.34900 1.97700 24.15000
H  30.78000 1.14300 24.68400
H  29.93300 1.13400 23.03700
H  30.14200 -0.30200 24.03500
H  31.60000 -2.42600 28.91000
H  32.38900 -0.93100 29.55000
H  31.71700 -1.01400 27.90200
H  29.35100 0.80400 33.71800
H  28.83300 -2.53300 36.59300
H  30.27100 -1.70300 37.18600
H  30.41000 -2.86500 36.01700
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


