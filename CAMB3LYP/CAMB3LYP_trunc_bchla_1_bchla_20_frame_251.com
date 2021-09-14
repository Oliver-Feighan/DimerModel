%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg -2.04400 17.49100 26.73800
C  -2.15500 15.53600 29.75200
C  -2.48400 20.32400 28.82400
C  -2.10000 19.56800 23.93100
C  -2.02900 14.73500 24.85000
N  -2.50300 17.81300 29.09300
C  -2.37100 16.88400 30.06400
C  -2.66300 17.53000 31.46600
C  -2.96600 19.01100 31.04000
C  -2.55400 19.10600 29.55300
C  -4.45800 19.41500 31.30700
C  -1.49000 17.42200 32.44000
C  -1.82900 17.03100 33.84100
H  -0.68700 16.80800 34.84300
N  -2.05000 19.60600 26.40500
C  -2.24700 20.56500 27.39300
C  -2.18900 21.87200 26.76900
C  -1.98200 21.71100 25.40900
C  -2.11500 20.23400 25.17200
C  -2.46700 23.16900 27.47300
C  -1.85900 22.77200 24.24000
O  -1.68000 22.55800 23.06600
C  -1.78900 24.24000 24.60100
N  -2.21300 17.23900 24.73600
C  -2.21800 18.17300 23.73300
C  -2.50700 17.54400 22.35000
C  -2.23600 16.00100 22.60900
C  -2.12000 15.99700 24.18100
C  -3.92800 17.89500 21.80400
C  -0.96100 15.44100 21.90800
C  0.37100 15.96600 22.43900
N  -2.11100 15.55600 27.12300
C  -2.02800 14.45800 26.22200
C  -2.06000 13.19200 26.95100
C  -2.05400 13.57600 28.30400
C  -2.07200 14.99400 28.37800
C  -2.02800 11.85900 26.36100
C  -2.05700 13.00500 29.64700
O  -1.98500 11.83800 29.99400
C  -2.21400 14.25700 30.66700
C  -1.01500 13.99300 31.58100
O  0.19500 14.13800 31.33000
O  -1.45400 13.50000 32.72600
C  -0.50700 13.41700 33.91100
H  -2.54500 21.19300 29.48200
H  -2.11900 20.13500 22.99700
H  -1.94800 13.88000 24.17600
H  -3.47200 16.92500 31.87400
H  -2.32200 19.72000 31.56200
H  -4.61000 20.30200 31.92200
H  -5.09400 18.58600 31.61800
H  -5.01600 19.78600 30.44700
H  -1.29200 18.47500 32.64100
H  -0.61800 16.92800 32.01000
H  -2.49700 16.18800 34.02500
H  -2.37000 17.90200 34.21100
H  -2.54100 23.02800 28.55100
H  -3.37300 23.68600 27.15300
H  -1.69300 23.93100 27.38000
H  -0.93800 24.47000 25.24200
H  -2.76500 24.39500 25.06100
H  -1.83600 24.73500 23.63100
H  -1.82000 18.06600 21.68300
H  -3.08800 15.32400 22.54300
H  -3.89400 18.55600 20.93800
H  -4.67200 18.26800 22.50700
H  -4.29000 16.90500 21.52600
H  -1.14700 15.45600 20.83400
H  -0.92500 14.37100 22.11200
H  0.30800 16.76300 23.17900
H  0.88600 16.37900 21.57200
H  0.96700 15.25500 23.01000
H  -2.06800 11.03700 27.07600
H  -1.09600 11.82400 25.79700
H  -2.85300 11.79600 25.65000
H  -3.07800 14.15100 31.32300
H  -1.14400 13.30000 34.78800
H  0.04800 14.34000 34.07600
H  0.23900 12.62700 33.81700
Mg 6.73600 57.36100 41.87100
C  5.57600 54.03200 41.84400
C  9.54700 56.26000 40.29200
C  7.50300 60.62700 41.21800
C  3.48900 58.25400 43.03700
N  7.42000 55.35200 41.01200
C  6.89900 54.10900 41.23100
C  7.68800 53.00900 40.50000
C  9.08300 53.80000 40.39600
C  8.66300 55.23700 40.58100
C  10.16800 53.39100 41.42800
C  7.12500 52.63700 39.05500
C  7.43400 51.26500 38.47600
H  8.38100 51.07500 37.27200
N  8.34400 58.41300 40.90600
C  9.44600 57.74300 40.46000
C  10.47100 58.62900 40.12400
C  9.93400 59.94200 40.40300
C  8.52300 59.74000 40.83800
C  11.84700 58.22300 39.61000
C  10.74100 61.20700 40.21700
O  11.93000 61.13100 39.81300
C  10.22200 62.63600 40.60600
N  5.63100 59.10700 42.11600
C  6.19700 60.37100 41.79700
C  5.36400 61.47100 42.39400
C  3.99000 60.81500 42.65100
C  4.37000 59.30600 42.68800
C  5.80500 62.13900 43.71600
C  2.94100 61.09300 41.45000
C  2.15100 62.39300 41.46200
N  4.97100 56.31000 42.51500
C  3.75500 56.80200 42.96500
C  2.84300 55.78700 43.18500
C  3.53400 54.62400 42.80400
C  4.81800 55.00700 42.38100
C  1.45400 55.97300 43.70900
C  3.46600 53.25900 42.78000
O  2.61500 52.50900 43.21600
C  4.68000 52.84700 41.92100
C  4.09400 52.40300 40.61400
O  3.41000 53.09000 39.92100
O  4.40200 51.06200 40.37800
C  3.73000 50.54700 39.23600
H  10.54100 55.93300 39.97900
H  7.77500 61.68500 41.22200
H  2.57200 58.59300 43.52300
H  7.89200 52.14900 41.13700
H  9.49500 53.56500 39.41500
H  10.79300 52.70800 40.85200
H  9.75000 52.84400 42.27300
H  10.76800 54.24600 41.74000
H  7.52400 53.32000 38.30500
H  6.05400 52.83900 39.08200
H  6.45300 50.86600 38.21400
H  7.81100 50.65000 39.29300
H  12.72600 58.79500 39.90900
H  11.62200 58.24700 38.54400
H  12.11100 57.23800 39.99500
H  10.03600 62.67400 41.67900
H  9.32800 62.73100 39.99000
H  10.83200 63.47600 40.27400
H  5.23800 62.21700 41.61000
H  3.54200 61.14600 43.58800
H  5.15700 61.94700 44.57100
H  6.23700 63.12400 43.53800
H  6.70300 61.58800 43.99500
H  2.28600 60.23100 41.57500
H  3.43700 61.03400 40.48100
H  2.26900 62.89100 40.50000
H  2.53500 63.03900 42.25200
H  1.08300 62.20800 41.58000
H  0.80400 55.41700 43.03400
H  1.12000 57.01000 43.65700
H  1.32300 55.57400 44.71500
H  5.17600 52.01400 42.41800
H  2.78900 50.09000 39.54100
H  4.23500 49.69300 38.78300
H  3.52200 51.31600 38.49400


