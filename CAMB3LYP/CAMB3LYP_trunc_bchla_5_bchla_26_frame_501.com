%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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


