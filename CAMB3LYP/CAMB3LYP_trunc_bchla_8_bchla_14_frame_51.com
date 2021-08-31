%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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
Mg 46.69300 44.58300 44.38700
C  43.25500 44.72400 43.60300
C  46.82600 41.39700 42.97300
C  50.08000 44.92600 44.22900
C  46.41600 48.03600 45.29800
N  45.21700 43.25300 43.45800
C  43.91500 43.52200 43.22400
C  43.24400 42.37600 42.39900
C  44.33300 41.22100 42.58000
C  45.59600 41.99100 43.08500
C  44.09700 40.05700 43.56900
C  43.05200 42.86300 40.91300
C  41.75200 42.34300 40.23900
H  41.89600 41.90400 38.72500
N  48.30200 43.23100 43.87300
C  48.09600 41.97500 43.36700
C  49.37200 41.38600 43.21500
C  50.37100 42.36600 43.49200
C  49.61800 43.58400 43.94300
C  49.62300 39.97500 42.73800
C  51.83400 42.25100 43.23000
O  52.31500 41.27900 42.67000
C  52.80700 43.26800 43.53300
N  48.03900 46.33200 44.62600
C  49.42400 46.11300 44.66400
C  50.13700 47.33600 45.15000
C  48.99300 48.37700 45.43900
C  47.74400 47.57600 45.11400
C  50.98200 46.96700 46.37600
C  49.23100 49.78300 44.77000
C  49.25400 50.97200 45.68200
N  45.14700 46.13000 44.49400
C  45.21100 47.43000 44.93200
C  43.85400 48.01300 45.01800
C  43.09500 46.96600 44.43100
C  43.89400 45.85700 44.22500
C  43.42700 49.41800 45.41300
C  41.78600 46.66700 43.91400
O  40.77400 47.38200 43.84700
C  41.78000 45.23600 43.37700
C  40.72900 44.42700 44.02200
O  39.86700 43.90600 43.41200
O  40.96200 44.18000 45.34400
C  40.02400 43.15900 45.85900
H  46.89000 40.36900 42.60900
H  51.17100 44.93400 44.24400
H  46.14300 49.07400 45.50200
H  42.28400 42.14500 42.86100
H  44.78800 40.86200 41.65700
H  43.41900 40.36700 44.36500
H  44.98600 39.61600 44.02000
H  43.58300 39.28600 42.99600
H  43.92400 42.61300 40.30800
H  43.01600 43.95000 40.98600
H  41.07000 43.18900 40.32600
H  41.42600 41.56000 40.92400
H  50.58400 39.53700 43.00800
H  49.53400 40.05600 41.65500
H  48.81600 39.25800 42.88900
H  52.77800 43.41800 44.61300
H  52.61000 44.18200 42.97200
H  53.79000 42.90400 43.23400
H  50.70700 47.79400 44.34200
H  48.71000 48.54100 46.47900
H  51.11500 45.89200 46.50000
H  50.52200 47.31200 47.30200
H  52.01900 47.25200 46.20100
H  48.42600 49.87600 44.04100
H  50.21400 49.89500 44.31400
H  48.51600 51.67900 45.30400
H  50.21500 51.48600 45.65500
H  49.14300 50.70900 46.73400
H  43.51500 50.13100 44.59300
H  43.93200 49.78600 46.30700
H  42.35600 49.27300 45.55400
H  41.46600 45.23800 42.33400
H  40.25400 43.00100 46.91300
H  40.12100 42.22900 45.29900
H  38.98500 43.48400 45.90100


