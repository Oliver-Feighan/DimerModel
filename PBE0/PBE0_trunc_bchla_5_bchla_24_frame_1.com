%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 24.45900 -7.80400 45.73000
C  26.58600 -5.04500 44.57700
C  21.88700 -6.07400 44.19000
C  22.84100 -10.55200 45.75100
C  27.52200 -9.37000 46.76800
N  24.34300 -5.76600 44.35600
C  25.30500 -4.83000 44.04500
C  24.68800 -3.74900 43.08400
C  23.16700 -3.94300 43.21100
C  23.11900 -5.33100 43.95600
C  22.33200 -2.87300 44.03000
C  25.23700 -3.82200 41.61800
C  25.48800 -2.49900 40.92700
H  24.41700 -2.03400 39.89100
N  22.60400 -8.23500 45.19100
C  21.64700 -7.40100 44.68000
C  20.38100 -8.15500 44.46900
C  20.67900 -9.47200 44.90600
C  22.09500 -9.47000 45.26500
C  19.08200 -7.63000 43.83600
C  19.81800 -10.73300 44.78900
O  18.67600 -10.59500 44.35900
C  20.22700 -12.02800 45.16700
N  25.20700 -9.77600 45.95500
C  24.19400 -10.69000 46.16700
C  24.84300 -12.02000 46.71400
C  26.38100 -11.70100 46.72900
C  26.35300 -10.16200 46.56900
C  24.30300 -12.52000 48.08200
C  27.32200 -12.42400 45.68300
C  27.76100 -13.86800 45.98200
N  26.63200 -7.37200 45.76700
C  27.69000 -8.05400 46.33400
C  28.89800 -7.26500 46.33900
C  28.50100 -6.05100 45.66300
C  27.12300 -6.13300 45.25800
C  30.21600 -7.74700 46.86200
C  28.99300 -4.84100 45.17500
O  30.13700 -4.41700 45.08200
C  27.84800 -4.18900 44.34700
C  27.67100 -2.81900 44.95900
O  27.62100 -1.78000 44.29600
O  27.39300 -2.80300 46.38500
C  26.88100 -1.51100 46.83900
H  21.04400 -5.51700 43.77700
H  22.37000 -11.50900 45.98500
H  28.34500 -9.94300 47.20000
H  24.83100 -2.76600 43.53500
H  22.69900 -4.12200 42.24300
H  21.79800 -2.30300 43.27000
H  22.98800 -2.17400 44.54800
H  21.67600 -3.27700 44.80100
H  24.55300 -4.27200 40.89800
H  26.19500 -4.34200 41.59800
H  26.49400 -2.54600 40.51000
H  25.50300 -1.72600 41.69600
H  18.22200 -8.02700 44.37600
H  19.15100 -7.95200 42.79700
H  19.01200 -6.54700 43.92800
H  20.35400 -12.03800 46.25000
H  21.09800 -12.40800 44.63300
H  19.40000 -12.62300 44.77800
H  24.71500 -12.87500 46.05100
H  26.80800 -11.85300 47.72000
H  25.17900 -12.58000 48.72800
H  23.87600 -13.51200 47.93800
H  23.57100 -11.76000 48.35800
H  28.20300 -11.78200 45.65400
H  26.82800 -12.52100 44.71600
H  27.32600 -14.00800 46.97200
H  28.83900 -13.92900 46.13800
H  27.42200 -14.64300 45.29500
H  30.43000 -8.65700 46.30300
H  30.17700 -8.08000 47.89900
H  31.04400 -7.08900 46.59800
H  28.22500 -4.06500 43.33200
H  26.34600 -1.90700 47.70300
H  26.35400 -0.97600 46.04900
H  27.73900 -0.96400 47.23000
Mg -0.18700 43.41100 24.93300
C  1.60700 43.07600 28.07200
C  -2.89800 42.36800 26.64400
C  -1.73200 43.65500 22.14800
C  2.88800 44.15900 23.34500
N  -0.58100 42.80100 27.23000
C  0.26100 42.84200 28.27700
C  -0.50600 42.58400 29.58800
C  -1.80900 41.92300 28.93700
C  -1.80000 42.43200 27.55700
C  -2.07500 40.44900 29.19300
C  -0.70200 43.92400 30.31800
C  -0.86900 43.75600 31.85700
H  -0.56100 42.38900 32.62600
N  -2.11500 43.10200 24.42900
C  -3.10300 42.58200 25.24300
C  -4.30800 42.41600 24.43600
C  -4.00300 42.86400 23.12700
C  -2.55900 43.22400 23.20800
C  -5.57100 41.87400 24.94000
C  -4.95600 43.00700 21.82800
O  -4.58400 43.57000 20.83500
C  -6.46400 42.70700 21.90000
N  0.55900 43.58500 22.97300
C  -0.33400 43.81600 22.01400
C  0.26300 44.22400 20.69400
C  1.62800 44.75300 21.19800
C  1.73900 43.97100 22.55200
C  0.44200 43.02300 19.65300
C  1.76200 46.32800 21.33100
C  0.94100 47.25900 20.39200
N  1.83500 43.48300 25.49800
C  2.95000 43.84600 24.81700
C  4.05200 43.97800 25.65800
C  3.57000 43.63900 26.89000
C  2.22100 43.36000 26.79400
C  5.44000 44.46500 25.25100
C  4.00800 43.52200 28.28900
O  5.10800 43.51800 28.82100
C  2.67500 43.22200 29.10400
C  2.92200 42.05100 29.92700
O  2.92900 40.87200 29.55900
O  3.40200 42.49600 31.10700
C  4.04200 41.44600 31.87900
H  -3.71800 41.86000 27.15600
H  -2.28400 43.68900 21.20700
H  3.79200 44.48700 22.82800
H  -0.02300 41.81500 30.19100
H  -2.65500 42.40200 29.42800
H  -2.95800 40.29300 29.81200
H  -1.21600 39.99100 29.68300
H  -2.18500 39.85300 28.28700
H  -1.65600 44.33500 29.98700
H  0.05400 44.67300 30.08300
H  -1.83500 44.12700 32.20100
H  -0.27300 44.58700 32.23400
H  -5.66200 41.99700 26.01900
H  -5.55700 40.82000 24.66200
H  -6.45800 42.34700 24.51800
H  -7.01900 43.39100 22.54200
H  -6.51800 41.64400 22.13600
H  -6.89500 42.73900 20.89900
H  -0.45600 44.95300 20.32100
H  2.46400 44.44400 20.57100
H  0.05400 43.47000 18.73800
H  -0.13700 42.13000 19.89100
H  1.46800 42.71200 19.46000
H  2.81700 46.59200 21.27300
H  1.31000 46.48200 22.31100
H  -0.04700 47.44600 20.81200
H  0.97500 46.73900 19.43400
H  1.44300 48.20900 20.21200
H  5.25700 45.13800 24.41300
H  6.07800 43.67100 24.86400
H  5.88500 44.97700 26.10400
H  2.48900 44.10200 29.72100
H  4.45500 41.79500 32.82500
H  4.72400 40.78900 31.33900
H  3.30000 40.72800 32.23000


