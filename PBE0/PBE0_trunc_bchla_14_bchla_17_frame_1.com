%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
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
Mg 29.46400 59.14500 41.03800
C  26.57800 57.38900 40.28900
C  31.34900 56.53100 39.78200
C  32.15000 61.18900 41.08600
C  27.35400 61.82000 41.94900
N  28.99200 57.29200 39.87100
C  27.77200 56.74800 39.89300
C  27.77900 55.38900 39.12700
C  29.32900 54.99700 39.21800
C  29.96700 56.34400 39.73200
C  29.58100 53.72700 40.06600
C  27.37400 55.59400 37.68500
C  26.46000 54.52700 37.10300
H  26.57800 54.02100 35.61600
N  31.53300 58.78800 40.65000
C  32.11000 57.71400 40.08400
C  33.51500 57.95800 40.03300
C  33.80600 59.24700 40.57400
C  32.47900 59.81600 40.84000
C  34.47700 56.89100 39.66100
C  35.16200 59.85800 40.59300
O  36.13000 59.21400 40.22400
C  35.39700 61.26400 41.13200
N  29.71300 61.28200 41.56500
C  30.94400 61.84600 41.41400
C  30.85900 63.33400 41.52000
C  29.29200 63.52200 41.79300
C  28.69900 62.13600 41.84900
C  31.79900 63.96200 42.56600
C  28.51800 64.48000 40.78300
C  27.64200 65.63900 41.31400
N  27.43700 59.54700 41.14000
C  26.74400 60.62700 41.58900
C  25.41200 60.28000 41.81900
C  25.25400 58.98400 41.34800
C  26.48900 58.61100 40.88000
C  24.32700 61.25500 42.36300
C  24.37200 57.92800 41.16800
O  23.15900 57.91800 41.34600
C  25.14900 56.75500 40.39700
C  25.25800 55.58200 41.18700
O  25.99300 55.31600 42.14100
O  24.33400 54.68200 40.66500
C  24.17500 53.35900 41.24400
H  31.95200 55.65500 39.53500
H  32.91900 61.96200 41.15300
H  26.80500 62.68200 42.33400
H  27.08300 54.71900 39.63100
H  29.81700 54.88800 38.25000
H  28.66900 53.37400 40.54700
H  30.34000 53.93800 40.82000
H  30.00600 52.97200 39.40500
H  28.16800 55.74800 36.95400
H  26.92700 56.58800 37.67000
H  25.44000 54.85100 37.30900
H  26.66400 53.69100 37.77200
H  34.81300 56.52500 40.63100
H  35.29700 57.28300 39.05900
H  34.02200 56.08300 39.08800
H  35.10600 61.93900 40.32700
H  36.47900 61.37300 41.20200
H  34.84400 61.42100 42.05800
H  31.14300 63.61000 40.50500
H  29.10900 64.04800 42.73000
H  32.59600 64.44300 41.99800
H  32.19100 63.23400 43.27700
H  31.26900 64.75100 43.09900
H  27.85100 63.86700 40.17700
H  29.18800 65.01800 40.11200
H  26.69000 65.11600 41.40400
H  27.54400 66.41700 40.55800
H  27.77900 66.00900 42.33000
H  24.06200 61.95200 41.56900
H  24.63400 61.90700 43.18100
H  23.42900 60.71200 42.65900
H  24.61600 56.61800 39.45600
H  23.13300 53.07400 41.10200
H  24.52900 53.25500 42.27000
H  24.72400 52.65900 40.61400


