%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 3.63500 0.41800 44.38500
C  6.38200 2.41600 44.15700
C  1.85600 2.91600 42.76500
C  1.03400 -1.70200 44.27100
C  5.65900 -2.27700 45.39800
N  4.15300 2.34100 43.32400
C  5.35700 2.93700 43.43600
C  5.44100 4.20100 42.66100
C  3.82600 4.53100 42.53700
C  3.19500 3.17000 42.89900
C  3.38200 5.72400 43.43400
C  6.08500 3.93300 41.31000
C  5.45900 2.81500 40.35400
H  5.36900 3.04800 38.85100
N  1.72200 0.57100 43.63100
C  1.11700 1.74400 43.11300
C  -0.27500 1.47100 42.97900
C  -0.55500 0.16000 43.34400
C  0.72000 -0.40500 43.81700
C  -1.17200 2.58100 42.47000
C  -1.87000 -0.57100 43.21500
O  -2.81800 0.09700 42.73200
C  -2.07300 -1.97000 43.39000
N  3.45400 -1.75300 44.57800
C  2.23500 -2.36400 44.54600
C  2.35900 -3.87900 44.83600
C  3.87500 -4.03200 45.19700
C  4.39100 -2.62000 44.97700
C  1.41200 -4.36800 45.91100
C  4.59000 -5.14900 44.26500
C  5.56000 -6.06500 45.02000
N  5.66900 0.04500 44.75000
C  6.32000 -0.99300 45.33600
C  7.63200 -0.52000 45.68500
C  7.75100 0.85800 45.29200
C  6.52700 1.13100 44.70400
C  8.77000 -1.38000 46.11700
C  8.57600 2.01500 45.15500
O  9.73800 2.10000 45.57200
C  7.71800 3.05600 44.46000
C  7.61900 4.24800 45.33700
O  8.63700 4.80700 45.69500
O  6.34700 4.66100 45.59200
C  6.23700 5.85100 46.43800
H  1.24500 3.78600 42.51600
H  0.13600 -2.29700 44.44900
H  6.31200 -3.00800 45.88000
H  6.00100 5.03200 43.09000
H  3.58600 4.93100 41.55100
H  2.93500 6.50500 42.82000
H  4.28700 6.13300 43.88400
H  2.70200 5.45700 44.24400
H  7.14700 3.76700 41.49200
H  6.12500 4.93200 40.87600
H  4.43400 2.52000 40.57900
H  6.17000 1.99900 40.48000
H  -1.10900 2.52700 41.38300
H  -0.76100 3.53900 42.78600
H  -2.21500 2.54500 42.78600
H  -2.88500 -2.34300 42.76500
H  -2.15000 -2.27300 44.43500
H  -1.20000 -2.43300 42.93200
H  2.14200 -4.33200 43.86800
H  3.97500 -4.21900 46.26600
H  0.70500 -3.61400 46.25600
H  1.77700 -4.91000 46.78300
H  0.73700 -5.01500 45.35000
H  5.15900 -4.73200 43.43400
H  3.86400 -5.76700 43.73700
H  5.09600 -7.04700 45.10800
H  5.65500 -5.69200 46.04000
H  6.52700 -6.16100 44.52600
H  8.71100 -1.83500 47.10500
H  9.65100 -0.74100 46.04600
H  8.91100 -2.23600 45.45700
H  8.25000 3.32500 43.54700
H  6.78900 6.65600 45.95200
H  6.55900 5.68400 47.46600
H  5.19100 6.09800 46.61900
Mg -9.42700 40.51200 41.96000
C  -8.27600 37.27100 41.01300
C  -7.18900 41.70700 39.82100
C  -11.06300 43.46900 41.99400
C  -12.24300 38.82400 43.34600
N  -7.88600 39.63200 40.71800
C  -7.50800 38.28500 40.52900
C  -6.30600 38.11700 39.54700
C  -5.96100 39.60100 39.23600
C  -7.01600 40.39000 40.04000
C  -4.48700 39.82100 39.61500
C  -6.54600 37.26600 38.22600
C  -5.77300 35.95800 38.10300
H  -4.88600 35.77200 36.85500
N  -9.14300 42.33400 41.14400
C  -8.11800 42.62800 40.24500
C  -8.09400 44.00500 39.95700
C  -9.10700 44.61100 40.67700
C  -9.83500 43.52300 41.29800
C  -7.10100 44.80200 39.03200
C  -9.33100 46.10500 40.68600
O  -8.59600 46.89800 40.17700
C  -10.57800 46.61100 41.56500
N  -11.40100 41.07400 42.51100
C  -11.80600 42.38000 42.46200
C  -13.34000 42.40200 42.85400
C  -13.69600 41.00800 43.39700
C  -12.37200 40.22400 43.05500
C  -13.75900 43.56700 43.70500
C  -14.84900 40.42700 42.58800
C  -15.68100 39.35600 43.24200
N  -10.03000 38.45700 42.37900
C  -11.16900 37.95500 42.96200
C  -11.08800 36.54300 43.06100
C  -9.91600 36.15600 42.36900
C  -9.38800 37.38300 41.87700
C  -12.09900 35.62000 43.61300
C  -9.03800 35.08300 41.88300
O  -9.11400 33.86600 42.09900
C  -7.91300 35.74600 40.99900
C  -6.65000 35.63000 41.67800
O  -6.37600 36.35300 42.66500
O  -5.96400 34.63500 41.08200
C  -4.53900 34.42700 41.45900
H  -6.43900 42.13800 39.15400
H  -11.53000 44.45500 42.04300
H  -13.01900 38.25900 43.86600
H  -5.46700 37.65700 40.06900
H  -6.10300 39.87100 38.19000
H  -4.30300 40.89300 39.55000
H  -3.75300 39.31500 38.98800
H  -4.26800 39.50000 40.63400
H  -6.11700 37.88500 37.43800
H  -7.61400 37.12100 38.06300
H  -6.50000 35.14700 38.15800
H  -5.08700 35.86800 38.94600
H  -7.73100 45.35500 38.33500
H  -6.44100 44.13800 38.47400
H  -6.37700 45.43500 39.54500
H  -10.52000 47.69900 41.55500
H  -10.39300 46.15100 42.53600
H  -11.42100 46.18000 41.02500
H  -13.83700 42.54600 41.89400
H  -13.95800 40.94300 44.45300
H  -14.59900 43.25300 44.32400
H  -14.03900 44.35800 43.00900
H  -12.95300 44.01500 44.28500
H  -14.53800 40.03800 41.61800
H  -15.54600 41.19600 42.25500
H  -15.47900 39.25300 44.30800
H  -15.45200 38.41100 42.74900
H  -16.73600 39.48600 43.00100
H  -12.69800 36.09600 44.39000
H  -11.53400 34.78400 44.02600
H  -12.74800 35.36800 42.77400
H  -7.89000 35.41100 39.96200
H  -4.47300 33.34300 41.35500
H  -4.28400 34.79600 42.45200
H  -3.86400 34.76500 40.67200


