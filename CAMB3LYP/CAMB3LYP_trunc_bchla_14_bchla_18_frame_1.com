%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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
Mg 35.49600 50.25400 25.27000
C  35.17700 48.23300 28.13300
C  33.95100 52.80200 27.06200
C  35.52400 52.00900 22.45100
C  36.25500 47.31600 23.35100
N  34.65000 50.47900 27.32000
C  34.75900 49.53200 28.34800
C  34.55800 50.22300 29.68400
C  33.90100 51.58900 29.33600
C  34.18900 51.63300 27.80500
C  32.41200 51.77300 29.76700
C  35.85100 50.33500 30.56400
C  35.57700 50.30300 32.13200
H  34.12200 50.22500 32.60200
N  34.96300 52.12400 24.87400
C  34.41800 53.04500 25.75200
C  34.29200 54.33000 25.12700
C  34.61900 54.10800 23.75100
C  35.08400 52.73500 23.66700
C  33.72300 55.59000 25.84800
C  34.52600 55.05800 22.60200
O  34.69900 54.67600 21.42000
C  34.14200 56.56700 22.60300
N  35.57000 49.71100 23.15300
C  35.66800 50.65100 22.19400
C  36.16600 50.02300 20.86800
C  36.37900 48.52900 21.13700
C  36.08600 48.52400 22.70200
C  35.30200 50.41100 19.70700
C  37.75400 47.93500 20.76300
C  38.97000 48.58500 21.46800
N  35.63500 48.14500 25.59300
C  36.09300 47.16800 24.77600
C  36.32500 45.99400 25.61600
C  36.01700 46.32200 26.96600
C  35.54300 47.65500 26.90400
C  36.93800 44.67500 25.15300
C  35.96700 45.83900 28.30900
O  36.34100 44.80600 28.85800
C  35.32300 47.11200 29.16700
C  34.10500 46.63100 29.93300
O  33.00700 46.49100 29.43100
O  34.39500 46.48000 31.27000
C  33.50200 45.67800 32.10000
H  33.57900 53.56200 27.75200
H  35.60400 52.66100 21.57900
H  36.61900 46.49900 22.72500
H  33.85700 49.64500 30.28600
H  34.49900 52.36000 29.82100
H  31.92700 52.02300 28.82400
H  32.30200 52.65000 30.40600
H  31.87100 50.93900 30.21400
H  36.34800 51.29400 30.42100
H  36.54700 49.55500 30.25800
H  35.89400 51.26500 32.53500
H  36.14900 49.47100 32.54200
H  33.11100 56.13200 25.12600
H  34.50500 56.27400 26.17600
H  33.12800 55.50300 26.75700
H  34.79000 57.12800 23.27700
H  33.09000 56.73000 22.83700
H  34.45500 57.11000 21.71200
H  37.07600 50.56600 20.61600
H  35.49800 48.04500 20.71400
H  35.74600 51.19400 19.09200
H  34.31200 50.69900 20.06100
H  35.07200 49.59000 19.02700
H  37.83100 48.14500 19.69600
H  37.75500 46.85700 20.92700
H  39.44400 47.77500 22.02200
H  38.57400 49.35600 22.12900
H  39.64500 49.11900 20.79900
H  37.67500 44.28600 25.85500
H  37.52900 44.80600 24.24700
H  36.12400 43.98800 24.92000
H  36.09300 47.35400 29.90000
H  32.48700 45.91500 31.78000
H  33.56400 45.88200 33.16900
H  33.81500 44.64500 31.94800


