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


