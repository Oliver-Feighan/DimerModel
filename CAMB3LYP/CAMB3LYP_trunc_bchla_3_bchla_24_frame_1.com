%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 1.39700 7.67200 26.33400
C  1.72400 9.74100 29.05100
C  2.01300 4.99100 28.38100
C  1.18200 5.57500 23.58000
C  1.22200 10.43500 24.24700
N  1.88100 7.36900 28.50400
C  1.86000 8.39700 29.40300
C  2.08800 7.91400 30.75300
C  2.40000 6.43300 30.60700
C  2.12200 6.20500 29.05700
C  3.81300 5.91300 31.04400
C  0.86900 8.21500 31.69000
C  1.18400 8.50600 33.19800
H  2.44800 7.85400 33.79200
N  1.46600 5.54900 26.07500
C  1.71900 4.62900 27.07800
C  1.70700 3.27200 26.49500
C  1.45400 3.42700 25.06400
C  1.33100 4.87900 24.84700
C  1.96100 1.99300 27.27600
C  1.22600 2.33500 24.04500
O  0.94200 2.65000 22.87600
C  1.34900 0.88300 24.32000
N  1.30400 7.94800 24.20100
C  1.17000 6.94200 23.33200
C  0.89900 7.47700 21.91600
C  0.64800 9.06900 22.12400
C  1.00600 9.14600 23.65000
C  2.06300 7.12000 20.94100
C  -0.80600 9.55500 21.87400
C  -1.00500 10.45600 20.64600
N  1.52400 9.70800 26.54100
C  1.45100 10.73400 25.59600
C  1.56500 12.05000 26.26300
C  1.69000 11.65900 27.61700
C  1.66900 10.24300 27.72900
C  1.49300 13.37100 25.67200
C  1.73800 12.19000 28.99500
O  1.68200 13.31300 29.42100
C  1.80500 10.95700 29.96500
C  2.97400 11.12000 30.84400
O  4.12200 10.98100 30.39200
O  2.68100 11.28600 32.16300
C  3.77100 11.25200 33.13300
H  2.21100 4.15900 29.06000
H  1.15400 4.90100 22.72100
H  1.03400 11.27000 23.56900
H  2.96600 8.42700 31.14700
H  1.62000 5.81700 31.05600
H  4.56300 6.57700 31.47500
H  4.27400 5.43500 30.18000
H  3.77100 5.07400 31.73900
H  0.12500 7.42500 31.58700
H  0.27300 9.01100 31.24400
H  0.33700 8.21300 33.81900
H  1.41300 9.55400 33.38900
H  2.02800 2.15000 28.35200
H  2.91300 1.54300 26.99300
H  1.08700 1.34400 27.22400
H  1.03600 0.37200 23.41000
H  0.74800 0.54700 25.16600
H  2.36800 0.65900 24.63600
H  -0.06900 7.08000 21.61000
H  1.35900 9.68000 21.56800
H  1.72500 6.67300 20.00600
H  2.74700 6.40500 21.39900
H  2.59200 8.06400 20.80800
H  -1.22600 10.10700 22.71400
H  -1.49800 8.72100 21.75500
H  -0.24900 10.22000 19.89700
H  -0.86300 11.48900 20.96200
H  -2.04500 10.34200 20.34300
H  1.26400 14.14800 26.40100
H  0.76200 13.44400 24.86700
H  2.50500 13.65900 25.38700
H  0.85200 11.05500 30.48600
H  4.47400 12.08100 33.04700
H  4.27000 10.28300 33.12600
H  3.36500 11.36600 34.13800
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


