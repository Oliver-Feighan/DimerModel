%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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


