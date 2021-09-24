%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 24.36900 -6.98700 46.60800
C  26.68400 -4.76900 45.35200
C  21.90100 -5.12400 45.15800
C  22.46500 -9.68400 46.72500
C  27.03400 -8.80400 47.99500
N  24.41600 -5.40800 45.00500
C  25.43000 -4.51800 44.85100
C  25.02900 -3.31400 44.05100
C  23.50400 -3.26400 44.39000
C  23.25600 -4.65400 44.95700
C  23.12700 -2.08400 45.33100
C  25.27600 -3.52100 42.47300
C  25.72300 -2.21200 41.67100
H  24.68300 -1.93800 40.61300
N  22.38400 -7.27700 46.19700
C  21.52700 -6.46900 45.53500
C  20.22900 -7.11300 45.37800
C  20.37200 -8.43000 45.86100
C  21.79800 -8.50600 46.28600
C  19.04300 -6.37400 44.71200
C  19.31200 -9.48200 45.68900
O  18.25700 -9.18500 45.19200
C  19.48300 -10.96700 46.06200
N  24.62900 -8.97800 47.39100
C  23.73100 -9.85500 47.30000
C  24.17900 -11.20600 47.72000
C  25.67400 -10.94700 48.15800
C  25.81200 -9.52100 47.84300
C  23.22200 -11.82700 48.79300
C  26.70500 -11.85900 47.34500
C  28.13900 -12.05100 47.95500
N  26.42100 -6.83200 46.76500
C  27.35000 -7.59400 47.45000
C  28.63900 -6.96000 47.42800
C  28.41400 -5.78200 46.63300
C  27.10800 -5.81900 46.18400
C  29.90100 -7.42400 48.08100
C  29.04200 -4.60800 46.13300
O  30.18200 -4.13000 46.19600
C  28.00800 -4.00400 45.12600
C  27.97900 -2.53300 45.37600
O  28.65400 -1.69700 44.75900
O  27.14200 -2.25800 46.42200
C  27.16400 -0.79200 46.68700
H  21.02900 -4.53900 44.85700
H  21.83100 -10.57200 46.69600
H  27.87200 -9.43200 48.30700
H  25.56100 -2.49300 44.53000
H  22.94300 -3.07800 43.47400
H  22.46400 -1.40000 44.80200
H  24.04300 -1.52500 45.52500
H  22.72200 -2.51300 46.24700
H  24.47600 -4.01300 41.92000
H  26.07600 -4.24200 42.30800
H  26.63900 -2.51400 41.16300
H  25.99300 -1.30800 42.21700
H  18.63800 -6.86100 43.82500
H  19.24900 -5.33100 44.46900
H  18.26300 -6.29000 45.46800
H  19.62500 -11.12900 47.13000
H  20.27800 -11.46900 45.51000
H  18.60200 -11.45500 45.64400
H  24.15900 -11.77400 46.79000
H  25.92100 -11.00100 49.21800
H  22.84600 -12.76100 48.37600
H  22.40400 -11.15300 49.04800
H  23.73300 -12.22400 49.67000
H  26.84400 -11.46000 46.34000
H  26.21500 -12.82700 47.25100
H  28.51200 -13.04100 48.21600
H  28.22300 -11.52500 48.90500
H  28.89100 -11.63100 47.28700
H  29.81500 -8.26400 48.77000
H  30.19100 -6.49700 48.57400
H  30.52400 -7.56700 47.19800
H  28.23100 -4.20700 44.07800
H  28.13300 -0.51500 47.10100
H  26.50700 -0.50200 47.50700
H  26.90400 -0.28600 45.75700
Mg 16.95800 -2.04300 28.05800
C  15.67900 0.03700 30.54100
C  18.78700 -3.72100 30.44000
C  18.31300 -3.78900 25.50500
C  15.32700 0.02400 25.64200
N  17.25900 -1.76300 30.19400
C  16.62800 -0.88400 31.01900
C  16.98400 -1.14100 32.48900
C  18.10000 -2.26900 32.35300
C  18.02400 -2.61400 30.86200
C  19.48200 -1.75800 32.74200
C  15.83300 -1.72700 33.38200
C  15.87500 -1.41800 34.85600
H  17.10800 -0.90300 35.45200
N  18.26300 -3.57700 28.00200
C  18.92400 -4.12800 29.11800
C  19.62000 -5.29300 28.63800
C  19.47300 -5.38800 27.23100
C  18.65200 -4.21500 26.84200
C  20.55500 -6.11300 29.50100
C  19.94500 -6.39800 26.27100
O  19.91000 -6.24200 25.06900
C  20.68300 -7.69100 26.77500
N  16.72600 -2.04100 25.82800
C  17.40600 -2.83200 25.05300
C  17.12300 -2.53800 23.54900
C  16.12900 -1.39100 23.56700
C  16.05500 -1.03300 25.09500
C  18.49500 -2.31000 22.75600
C  14.77700 -1.75000 22.91500
C  14.57500 -1.34500 21.40700
N  15.59300 -0.54700 28.03400
C  15.01200 0.20800 26.98500
C  14.13500 1.19500 27.48900
C  14.45500 1.19400 28.91300
C  15.31700 0.08300 29.13900
C  13.41400 2.25000 26.67000
C  14.16900 1.92500 30.09200
O  13.43900 2.88400 30.30900
C  14.95100 1.07800 31.29700
C  15.78000 2.10000 32.01800
O  16.65500 2.79200 31.48700
O  15.36200 2.29400 33.36800
C  16.13600 3.30900 34.14900
H  19.34000 -4.20600 31.24800
H  18.93000 -4.26900 24.74300
H  14.76600 0.65500 24.95000
H  17.39400 -0.20500 32.86900
H  17.90900 -3.13400 32.98800
H  19.41600 -0.81500 33.28400
H  20.16500 -1.70100 31.89400
H  19.97000 -2.35300 33.51400
H  15.82100 -2.80600 33.23200
H  14.89600 -1.40800 32.92500
H  15.72700 -2.36200 35.38000
H  15.08500 -0.70500 35.09200
H  21.59900 -6.03400 29.19800
H  20.20800 -7.14700 29.48400
H  20.60300 -5.83600 30.55400
H  21.71600 -7.39100 26.95300
H  20.76400 -8.41600 25.96500
H  20.25400 -8.04400 27.71200
H  16.67800 -3.47000 23.20000
H  16.55100 -0.48300 23.13900
H  18.49400 -1.34000 22.25900
H  18.63800 -3.00100 21.92500
H  19.26300 -2.34000 23.53000
H  14.01500 -1.24300 23.50800
H  14.64200 -2.82300 23.05200
H  14.22000 -0.31400 21.38800
H  13.95900 -1.97800 20.76800
H  15.61200 -1.30800 21.07500
H  12.51000 1.73500 26.34600
H  13.99600 2.54300 25.79700
H  13.31000 3.23400 27.12900
H  14.20300 0.59900 31.93000
H  16.03100 4.30700 33.72500
H  17.20500 3.14900 34.29200
H  15.57400 3.26700 35.08200


