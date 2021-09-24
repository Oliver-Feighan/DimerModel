%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg 40.54200 8.71200 29.57000
C  42.46900 10.35300 31.83700
C  38.24000 7.70200 31.93300
C  38.87500 6.96300 27.14400
C  43.10200 9.33000 27.08900
N  40.50700 8.83000 31.66000
C  41.28100 9.79100 32.39300
C  40.70700 10.01900 33.73400
C  39.62800 8.84800 33.81600
C  39.41100 8.39300 32.39000
C  40.09700 7.61800 34.74300
C  40.24300 11.48300 34.00600
C  40.64100 12.04200 35.33400
H  40.55300 11.05000 36.56000
N  38.75400 7.55700 29.53000
C  38.00600 7.28200 30.63600
C  36.89200 6.50300 30.21600
C  36.97400 6.32300 28.78000
C  38.23000 7.03300 28.36500
C  35.87900 5.90900 31.21900
C  36.08100 5.65900 27.77900
O  36.28700 5.46300 26.57600
C  34.68000 5.19900 28.33400
N  40.98100 8.10600 27.39000
C  40.09800 7.49500 26.64000
C  40.46000 7.54300 25.15700
C  41.81600 8.30500 25.16100
C  42.00400 8.59400 26.64900
C  40.56900 6.08800 24.47000
C  41.66000 9.63200 24.29500
C  40.74000 10.78700 24.76200
N  42.48900 9.50900 29.49800
C  43.34000 9.72700 28.43100
C  44.44800 10.55700 28.94600
C  44.20300 10.86700 30.25000
C  42.96700 10.18300 30.54000
C  45.69600 10.84800 28.18600
C  44.66700 11.56800 31.44700
O  45.66700 12.23600 31.65000
C  43.53700 11.23400 32.51700
C  44.17900 10.46300 33.60100
O  44.82400 9.40400 33.50900
O  43.90700 11.19500 34.73900
C  44.28800 10.53200 35.99300
H  37.46200 7.40700 32.64000
H  38.23300 6.48700 26.39900
H  43.84700 9.59300 26.33600
H  41.53500 9.82500 34.41700
H  38.70200 9.27400 34.20100
H  40.95100 7.85100 35.37900
H  40.49700 6.86600 34.06300
H  39.17400 7.30300 35.22800
H  39.15600 11.54000 33.94800
H  40.73700 12.04000 33.21000
H  39.99600 12.85500 35.66500
H  41.60800 12.52200 35.48500
H  35.14900 6.60800 31.62700
H  36.53200 5.46300 31.96900
H  35.20300 5.18600 30.76300
H  34.82100 4.32000 28.96300
H  33.98600 5.15200 27.49500
H  34.17100 5.96300 28.92200
H  39.71300 8.14400 24.63900
H  42.64100 7.73100 24.73900
H  39.67900 5.98600 23.84800
H  40.55600 5.31000 25.23400
H  41.39700 5.86900 23.79700
H  41.39900 9.47200 23.24900
H  42.67000 10.03800 24.36000
H  41.45000 11.59500 24.93900
H  40.14900 10.52600 25.63900
H  40.01100 11.10100 24.01500
H  45.97000 11.83700 28.55200
H  45.54800 10.83600 27.10600
H  46.48500 10.23100 28.61400
H  43.32700 12.24400 32.86800
H  44.09900 11.31400 36.72900
H  45.36500 10.37600 35.93100
H  43.84100 9.55300 36.16700


