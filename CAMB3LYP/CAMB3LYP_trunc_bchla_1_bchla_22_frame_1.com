%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg -1.96800 16.99100 27.07800
C  -2.18500 14.77300 29.87900
C  -3.22300 19.43300 29.10400
C  -2.21300 18.84100 24.51400
C  -1.82200 14.16700 25.04100
N  -2.58500 16.98200 29.20400
C  -2.56200 16.08300 30.18300
C  -2.90400 16.68200 31.58000
C  -3.74100 17.95500 31.11000
C  -3.12300 18.15900 29.73800
C  -5.23100 17.67600 30.94100
C  -1.62300 17.10100 32.37700
C  -1.84400 17.38800 33.86500
H  -0.80000 17.02400 34.81100
N  -2.15300 18.95200 26.92600
C  -2.67100 19.79100 27.85500
C  -2.58800 21.14700 27.27100
C  -2.15500 21.00600 25.85300
C  -2.18200 19.57400 25.69800
C  -2.85500 22.39700 28.05200
C  -1.78700 22.00400 24.80100
O  -1.66500 21.68100 23.57000
C  -1.65800 23.40500 25.24100
N  -2.23300 16.57800 25.04900
C  -2.24700 17.56000 24.14100
C  -2.16400 17.01500 22.74600
C  -1.84400 15.45800 22.89800
C  -1.88200 15.35900 24.44300
C  -3.41900 17.19600 21.88100
C  -0.59500 14.92100 22.15500
C  0.65800 15.50100 22.77200
N  -1.81800 14.84700 27.38100
C  -1.79500 13.85500 26.39000
C  -1.79500 12.55300 27.08200
C  -1.83500 12.88000 28.43600
C  -1.88500 14.29600 28.58000
C  -1.81900 11.27600 26.45500
C  -1.87200 12.36200 29.79400
O  -1.76000 11.24400 30.29900
C  -2.22300 13.56300 30.77900
C  -1.25100 13.50900 31.92600
O  -0.03700 13.73600 31.74700
O  -1.93000 13.43500 33.08100
C  -1.02500 13.40800 34.26100
H  -3.82700 20.21100 29.57500
H  -2.18600 19.44400 23.60400
H  -1.69200 13.24400 24.47300
H  -3.52200 16.04200 32.21000
H  -3.55300 18.82200 31.74300
H  -5.45200 16.61400 31.04700
H  -5.47900 18.01800 29.93600
H  -5.82000 18.29200 31.62100
H  -1.29900 18.08700 32.04400
H  -0.83800 16.35200 32.26700
H  -2.81200 17.01700 34.20300
H  -1.90600 18.46800 33.99900
H  -3.31200 22.17700 29.01600
H  -3.61900 22.95400 27.51000
H  -1.96100 23.01100 28.16200
H  -2.60300 23.73500 25.67300
H  -1.30300 23.96100 24.37400
H  -0.92000 23.45700 26.04200
H  -1.36000 17.49200 22.18700
H  -2.66000 14.83500 22.53100
H  -3.23500 17.73100 20.94900
H  -4.17200 17.70500 22.48200
H  -3.93500 16.27100 21.62500
H  -0.61200 15.23200 21.11100
H  -0.51600 13.83400 22.17600
H  0.74600 16.46900 22.28000
H  1.51900 14.87000 22.54900
H  0.52700 15.68800 23.83800
H  -1.10100 10.65800 26.99400
H  -1.56400 11.34900 25.39800
H  -2.85900 10.99200 26.61500
H  -3.25100 13.41000 31.10600
H  -1.27000 12.51300 34.83300
H  -1.13500 14.37200 34.75800
H  0.03200 13.23700 34.05400
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


