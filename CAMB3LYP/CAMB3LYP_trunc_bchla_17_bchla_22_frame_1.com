%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
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


