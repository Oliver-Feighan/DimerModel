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


