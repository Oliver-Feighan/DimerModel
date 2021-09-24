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
Mg -2.33100 34.27800 27.06700
C  -3.29300 32.73500 29.99200
C  -0.77900 36.82100 29.04000
C  -1.68900 35.99300 24.34500
C  -4.00500 31.86100 25.19100
N  -2.10600 34.64700 29.26400
C  -2.54800 33.87100 30.25800
C  -2.04400 34.45500 31.61800
C  -1.35700 35.84600 31.23800
C  -1.46000 35.80800 29.70100
C  -2.07400 37.11000 31.67200
C  -1.16200 33.50400 32.44100
C  -1.65900 33.38600 33.87200
H  -0.69900 33.03500 35.01800
N  -1.29600 36.10900 26.75300
C  -0.70200 36.96400 27.63600
C  -0.06500 38.02800 26.99800
C  -0.40300 37.90700 25.60500
C  -1.15900 36.61100 25.50300
C  0.74500 39.09500 27.72000
C  -0.12100 38.87200 24.46500
O  -0.65300 38.65500 23.35900
C  0.65800 40.10800 24.66200
N  -2.76200 33.93100 25.11000
C  -2.33900 34.76200 24.09300
C  -2.78500 34.31300 22.65200
C  -3.27700 32.84000 22.98500
C  -3.37800 32.87200 24.52100
C  -3.86900 35.12000 21.83400
C  -2.50400 31.61000 22.32400
C  -1.18700 31.25300 22.96100
N  -3.52800 32.72900 27.43000
C  -4.16000 31.81100 26.61300
C  -4.85300 30.83200 27.34600
C  -4.49800 31.15000 28.66200
C  -3.68000 32.31700 28.68800
C  -5.57400 29.72900 26.69600
C  -4.80700 30.73800 29.97700
O  -5.59800 29.85900 30.32500
C  -4.00900 31.74000 30.96400
C  -4.95800 32.38400 31.85500
O  -5.73200 33.27000 31.57300
O  -4.95600 31.78800 33.04100
C  -5.91200 32.19500 34.02600
H  -0.35300 37.61900 29.65100
H  -1.55300 36.44800 23.36200
H  -4.51800 31.12700 24.56600
H  -2.99900 34.65700 32.10200
H  -0.41300 35.85600 31.78400
H  -2.94800 36.94200 32.30200
H  -2.45800 37.64200 30.80200
H  -1.35600 37.81800 32.08700
H  -0.10200 33.75800 32.45800
H  -1.18900 32.52900 31.95500
H  -2.34600 32.53900 33.87800
H  -2.14000 34.33700 34.10000
H  0.94800 38.81600 28.75400
H  0.13200 39.99100 27.82000
H  1.62000 39.38400 27.13900
H  0.44000 40.83300 23.87700
H  1.70900 39.83900 24.55200
H  0.43200 40.54800 25.63300
H  -1.93500 34.12600 21.99500
H  -4.26000 32.69700 22.53800
H  -4.17100 35.97400 22.43900
H  -4.79700 34.54700 21.82100
H  -3.62200 35.43900 20.82100
H  -2.25200 31.92600 21.31200
H  -3.11500 30.71700 22.19800
H  -1.06000 30.17800 22.83800
H  -0.96900 31.52500 23.99300
H  -0.42200 31.79500 22.40500
H  -5.20900 29.35700 25.73900
H  -6.66000 29.81200 26.67700
H  -5.43800 28.93700 27.43200
H  -3.30800 31.14200 31.54600
H  -5.54400 33.01100 34.64900
H  -6.05600 31.34600 34.69400
H  -6.83100 32.52000 33.53800


