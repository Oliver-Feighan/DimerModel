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
Mg -9.66500 18.33400 42.78400
C  -6.24000 17.80500 42.49600
C  -9.09900 21.67600 42.05300
C  -12.91200 18.87400 42.48300
C  -10.14800 14.91400 43.00500
N  -7.91600 19.56100 42.17100
C  -6.59800 19.16200 42.25800
C  -5.63300 20.29600 42.07600
C  -6.57300 21.55800 42.31900
C  -7.95100 20.93600 42.13600
C  -6.37100 22.54100 43.56800
C  -4.81400 20.28500 40.70300
C  -5.51800 19.69700 39.44400
H  -5.27900 20.48100 38.17600
N  -10.85700 20.01500 42.21400
C  -10.42000 21.32200 41.96500
C  -11.51000 22.18000 41.58900
C  -12.66600 21.36900 41.68400
C  -12.20600 20.01400 42.11100
C  -11.29500 23.68900 41.46200
C  -14.07900 21.93000 41.40500
O  -14.25800 23.08100 40.99400
C  -15.20100 21.03000 41.71200
N  -11.24800 17.01300 42.63400
C  -12.52300 17.52900 42.79400
C  -13.59000 16.43700 43.02800
C  -12.64200 15.17100 42.94600
C  -11.28500 15.70200 42.80600
C  -14.41600 16.47100 44.27500
C  -12.97300 14.21900 41.68800
C  -13.02000 12.77100 42.10100
N  -8.45700 16.64500 42.95800
C  -8.80300 15.30900 43.04100
C  -7.58100 14.47400 43.08900
C  -6.53000 15.43300 42.87800
C  -7.12900 16.69600 42.80500
C  -7.42000 13.00900 43.03700
C  -5.08700 15.65900 42.64900
O  -4.19400 14.85200 42.53200
C  -4.84400 17.18800 42.43700
C  -3.86600 17.68600 43.40500
O  -2.75700 18.18900 43.22300
O  -4.28400 17.22800 44.68200
C  -3.45300 17.67600 45.82000
H  -8.86100 22.74200 42.07500
H  -13.99400 18.96100 42.60100
H  -10.40600 13.86400 43.15400
H  -4.85200 20.25200 42.83500
H  -6.35600 22.13700 41.42100
H  -5.49500 22.34800 44.18700
H  -7.18200 22.50000 44.29500
H  -6.14100 23.54200 43.20200
H  -3.84000 19.81700 40.84200
H  -4.55300 21.32400 40.50200
H  -6.59000 19.59400 39.61500
H  -5.14100 18.69400 39.24600
H  -11.94900 24.16700 40.73200
H  -10.29700 23.87500 41.06600
H  -11.49900 24.10000 42.45100
H  -15.27200 20.92300 42.79400
H  -15.05600 20.07500 41.20600
H  -16.08400 21.54300 41.33100
H  -14.18800 16.36800 42.11900
H  -12.69300 14.77700 43.96100
H  -14.92900 15.51000 44.31600
H  -15.20300 17.20800 44.11500
H  -13.79200 16.69500 45.14100
H  -12.21000 14.42700 40.93800
H  -13.86700 14.45400 41.10900
H  -12.15800 12.27100 41.66000
H  -13.81800 12.22200 41.60000
H  -12.89800 12.69100 43.18100
H  -8.34400 12.43000 43.00100
H  -6.93100 12.73100 43.97000
H  -6.96900 12.89900 42.05100
H  -4.46800 17.26400 41.41700
H  -4.21000 17.98200 46.54300
H  -2.74000 18.49300 45.70900
H  -3.06600 16.75200 46.25000


