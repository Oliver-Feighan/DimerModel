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
Mg 9.20000 47.36600 24.88000
C  7.02600 47.76900 27.74700
C  11.57000 48.44700 26.79100
C  10.94600 47.19300 22.09000
C  6.29100 46.96700 22.73800
N  9.27600 47.92100 27.10800
C  8.32000 48.03400 28.09600
C  8.93500 48.52700 29.35600
C  10.25500 49.18700 28.88100
C  10.38000 48.52600 27.48800
C  10.19000 50.72500 28.75600
C  9.27100 47.32000 30.28400
C  8.91000 47.48600 31.81100
H  8.76300 46.21000 32.73200
N  11.12300 47.75000 24.52700
C  11.99300 48.08000 25.44400
C  13.30700 48.11300 24.83000
C  13.09300 47.62600 23.51900
C  11.67000 47.53200 23.33600
C  14.56100 48.44000 25.66100
C  14.15200 47.33300 22.49300
O  13.84400 46.95300 21.40300
C  15.59100 47.42400 22.88400
N  8.59400 47.38800 22.62900
C  9.59400 47.15300 21.78200
C  9.10700 46.72400 20.44900
C  7.56400 46.85800 20.64800
C  7.41500 47.09900 22.12700
C  9.70100 47.44000 19.18300
C  6.80800 45.57800 20.11700
C  6.00900 45.65200 18.78600
N  7.17400 47.35700 25.08600
C  6.10400 47.11400 24.19600
C  4.83100 47.11900 24.94600
C  5.17200 47.34400 26.27300
C  6.54300 47.52700 26.33600
C  3.42100 47.05000 24.31100
C  4.55900 47.48400 27.63300
O  3.43200 47.38700 28.05800
C  5.73500 47.94200 28.63300
C  5.68600 47.23700 29.86100
O  6.32700 46.19300 30.04600
O  4.84600 47.74900 30.84500
C  4.90400 47.16700 32.22100
H  12.46500 48.75400 27.33500
H  11.53100 47.07200 21.17600
H  5.40800 46.83900 22.11000
H  8.32300 49.32300 29.78100
H  11.13800 49.00400 29.49300
H  10.92100 51.19300 29.41600
H  9.18200 51.02300 29.04300
H  10.33700 51.10600 27.74500
H  10.35100 47.19900 30.19800
H  8.78700 46.45500 29.83100
H  8.00100 48.08400 31.87000
H  9.74500 48.08400 32.17800
H  15.18900 47.54900 25.70000
H  14.38200 48.83100 26.66200
H  15.20000 49.13500 25.11700
H  15.85600 46.71400 23.66700
H  15.91500 48.35300 23.35300
H  16.20400 47.15000 22.02500
H  9.34100 45.66500 20.34000
H  7.22100 47.78300 20.18400
H  10.49600 48.09000 19.54700
H  9.00000 47.94700 18.51900
H  10.12400 46.67900 18.52700
H  6.07000 45.27200 20.85900
H  7.53000 44.76900 20.22800
H  6.23700 44.78300 18.16900
H  6.25500 46.57800 18.26600
H  4.93600 45.62800 18.98000
H  3.49200 46.43400 23.41500
H  3.14200 48.04200 23.95500
H  2.68000 46.64400 25.00000
H  5.61100 49.00700 28.82700
H  5.88500 47.32300 32.67200
H  4.47300 46.16700 32.25100
H  4.32200 47.83100 32.86000


