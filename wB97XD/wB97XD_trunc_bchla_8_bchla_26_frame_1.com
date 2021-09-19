%nproc=24
%mem=175gb
#p wB97XD/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 44.80500 2.87200 46.67800
C  42.56400 5.57100 46.81800
C  42.13800 0.79200 46.35700
C  47.01300 0.52400 46.12300
C  47.34200 5.37800 46.17800
N  42.64100 3.17100 46.66900
C  41.94200 4.34400 46.71700
C  40.48300 4.04900 46.33800
C  40.34000 2.59000 46.75700
C  41.78900 2.14300 46.48400
C  39.96100 2.42600 48.24700
C  40.23600 4.25700 44.79800
C  39.26900 5.45400 44.33400
H  39.83900 6.57400 43.44700
N  44.60700 0.93700 46.23600
C  43.41900 0.24600 46.16900
C  43.72000 -1.14400 46.02500
C  45.14200 -1.32500 46.06900
C  45.67600 0.06600 46.14000
C  42.62100 -2.25300 45.89000
C  45.87200 -2.71700 46.06400
O  45.23700 -3.72500 46.11700
C  47.31400 -2.74100 46.27700
N  46.83800 2.95400 45.90000
C  47.54500 1.82700 46.15300
C  49.01900 2.09000 46.41200
C  49.17500 3.56200 46.06200
C  47.67600 4.02200 46.03300
C  49.44600 1.78300 47.89500
C  49.96600 3.97900 44.73400
C  51.33300 4.60700 44.78400
N  45.00300 4.99200 46.72600
C  46.07700 5.86500 46.51000
C  45.65700 7.24900 46.55200
C  44.29400 7.16900 46.69900
C  43.89900 5.79300 46.77900
C  46.53700 8.41600 46.37700
C  43.08200 7.96600 46.81100
O  42.95000 9.19400 46.88400
C  41.94500 6.98400 47.00000
C  41.32800 7.26400 48.36600
O  42.00600 7.20100 49.38900
O  40.00100 7.70700 48.25300
C  39.36700 8.02300 49.50500
H  41.26200 0.14300 46.41400
H  47.71200 -0.30200 46.26300
H  48.13300 6.13000 46.13700
H  39.87300 4.80700 46.83000
H  39.59400 2.12600 46.11300
H  39.12500 1.73800 48.12600
H  39.66800 3.42600 48.56800
H  40.73100 2.03300 48.91100
H  39.76200 3.35000 44.42300
H  41.21000 4.30800 44.31200
H  38.98200 5.96200 45.25500
H  38.33700 5.13200 43.87100
H  42.71600 -2.97100 45.07600
H  41.68000 -1.71700 45.76200
H  42.51800 -2.85500 46.79300
H  47.88700 -2.17400 45.54300
H  47.64600 -3.77900 46.25800
H  47.54400 -2.49600 47.31300
H  49.67000 1.44400 45.82200
H  49.62700 4.10200 46.89400
H  49.91800 0.80800 48.01500
H  48.60400 1.96800 48.56200
H  50.20200 2.52100 48.16700
H  49.41200 4.48800 43.94600
H  50.09500 2.95900 44.37100
H  52.11800 3.97300 44.37300
H  51.60000 4.81700 45.82000
H  51.22600 5.52900 44.21300
H  47.39900 8.09600 45.79100
H  46.88900 8.66700 47.37800
H  45.99500 9.26600 45.96500
H  41.25100 7.27600 46.21300
H  39.04600 7.15600 50.08300
H  38.48400 8.64500 49.36100
H  39.98200 8.61300 50.18400
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


