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
Mg 53.03600 23.47300 43.74200
C  50.27300 25.68400 43.34900
C  50.89700 20.80500 43.54400
C  55.76900 21.39500 43.06400
C  55.19900 26.21400 43.34600
N  50.76700 23.35900 43.20300
C  49.85000 24.33200 43.33000
C  48.49000 23.71000 43.58900
C  48.66000 22.11500 43.58600
C  50.27200 22.11900 43.58300
C  48.04000 21.54600 44.88200
C  47.47500 24.17900 42.55900
C  47.98000 24.34500 41.10600
H  47.51400 23.28000 40.07600
N  53.32100 21.33300 43.43700
C  52.26200 20.45300 43.38400
C  52.84400 19.14800 43.34900
C  54.26800 19.24400 43.11100
C  54.55100 20.69600 43.23900
C  52.01300 17.88700 43.29400
C  55.22500 18.13300 42.67600
O  54.80000 17.01500 42.45700
C  56.63100 18.33100 42.48800
N  55.16000 23.75300 43.17000
C  56.04500 22.76900 43.10900
C  57.45600 23.24500 42.80400
C  57.16400 24.78400 42.53600
C  55.80000 24.97000 43.09600
C  58.49900 22.99600 44.00000
C  57.34100 25.12500 40.98500
C  58.47200 25.99800 40.52000
N  52.84200 25.47200 43.66600
C  53.82600 26.48100 43.63200
C  53.20700 27.81500 43.74300
C  51.85100 27.50000 43.65700
C  51.62900 26.08900 43.48500
C  53.87500 29.20800 44.01000
C  50.52800 28.13200 43.63200
O  50.14800 29.31100 43.72400
C  49.45600 27.00700 43.45200
C  48.88400 27.41400 42.13800
O  49.50700 27.53100 41.11400
O  47.60700 27.75600 42.28300
C  46.79800 28.12100 41.06800
H  50.22500 19.97800 43.78200
H  56.71400 20.86800 42.91800
H  55.94300 26.97900 43.11200
H  48.29500 24.16100 44.56200
H  48.23100 21.65000 42.69800
H  47.24600 20.85400 44.60200
H  47.57500 22.29200 45.52700
H  48.83900 20.99200 45.37400
H  47.06700 25.14300 42.86300
H  46.65900 23.46100 42.47200
H  49.05000 24.24200 40.92600
H  47.63400 25.33600 40.81200
H  50.94600 18.10600 43.24600
H  52.31600 17.25700 44.13100
H  52.35200 17.31500 42.43100
H  56.90400 19.01600 41.68500
H  57.10100 17.37500 42.25800
H  56.94000 18.76400 43.44000
H  57.86000 22.86500 41.86600
H  58.02800 25.24600 43.01400
H  58.00400 22.45500 44.80700
H  58.78500 23.97500 44.38500
H  59.43300 22.50300 43.72600
H  56.50400 25.71600 40.61300
H  57.40500 24.21500 40.39000
H  58.17800 27.02600 40.30800
H  58.72900 25.47900 39.59600
H  59.39500 25.89600 41.09100
H  54.34600 29.53100 43.08200
H  54.56000 29.07800 44.84800
H  53.08700 29.93100 44.22000
H  48.79600 27.08200 44.31700
H  46.90900 29.16800 40.78500
H  45.73600 27.93300 41.23100
H  47.15900 27.45300 40.28600


