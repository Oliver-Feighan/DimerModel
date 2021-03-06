%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg -2.17200 16.97900 26.95200
C  -2.72400 15.15000 29.85500
C  -3.03700 19.79700 28.58300
C  -2.39300 18.59300 23.98800
C  -2.09700 13.93400 25.14200
N  -2.70100 17.33200 28.98300
C  -2.73500 16.51500 30.09200
C  -3.14200 17.32600 31.35300
C  -3.73400 18.61600 30.70600
C  -3.02700 18.64900 29.35900
C  -5.29100 18.59800 30.51100
C  -1.92300 17.38000 32.29200
C  -2.24700 17.17600 33.72300
H  -1.11900 16.77000 34.55900
N  -2.42900 18.95600 26.42700
C  -2.71700 19.97700 27.27500
C  -2.49100 21.22100 26.58300
C  -2.29300 20.85800 25.20200
C  -2.45300 19.40200 25.14600
C  -2.60700 22.55600 27.22700
C  -2.12700 21.73600 23.91900
O  -1.95600 21.25500 22.79700
C  -2.02300 23.23300 24.00700
N  -2.58000 16.29000 24.86100
C  -2.48200 17.22200 23.79200
C  -2.41100 16.57800 22.44700
C  -2.09800 15.11800 22.91000
C  -2.20000 15.11800 24.43300
C  -3.71000 16.74000 21.58700
C  -0.67500 14.61400 22.47000
C  0.50900 15.54800 22.82800
N  -2.25800 14.94800 27.37500
C  -2.16300 13.82100 26.54400
C  -2.25500 12.68900 27.38800
C  -2.44400 13.17700 28.67500
C  -2.44900 14.53100 28.60300
C  -2.11400 11.26600 26.95000
C  -2.70100 12.78200 29.98000
O  -2.81800 11.63600 30.45100
C  -2.89700 14.02900 30.89400
C  -1.76000 13.88300 31.95300
O  -0.54300 14.00600 31.72700
O  -2.34900 13.58900 33.17800
C  -1.46500 13.15200 34.30500
H  -3.42500 20.61700 29.19100
H  -2.30900 18.94200 22.95600
H  -1.87100 12.97600 24.66900
H  -3.88000 16.77500 31.93500
H  -3.40400 19.48400 31.27600
H  -5.81600 17.73100 30.91300
H  -5.54300 18.53800 29.45300
H  -5.71700 19.45700 31.02900
H  -1.55900 18.40500 32.23700
H  -1.19600 16.65200 31.93100
H  -2.95800 16.35000 33.74700
H  -2.76300 18.06300 34.09000
H  -2.53500 22.48300 28.31200
H  -3.46700 23.06800 26.79400
H  -1.74500 23.16300 26.95000
H  -2.35600 23.58400 23.03000
H  -0.98600 23.54300 24.13800
H  -2.68700 23.64600 24.76700
H  -1.57800 16.92500 21.83400
H  -2.85900 14.45200 22.50200
H  -4.54500 16.97100 22.24900
H  -4.01900 15.87700 20.99700
H  -3.69300 17.53600 20.84200
H  -0.75400 14.40000 21.40400
H  -0.59800 13.63700 22.94700
H  1.24700 15.06400 23.46800
H  0.07900 16.42200 23.31800
H  1.05200 15.89400 21.94900
H  -3.02600 10.76600 27.27600
H  -1.20700 10.86800 27.40500
H  -2.07500 11.10900 25.87200
H  -3.86000 13.93400 31.39500
H  -0.97000 12.24300 33.96400
H  -1.93600 12.93200 35.26300
H  -0.82600 14.03500 34.31800
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


