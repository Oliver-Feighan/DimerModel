%nproc=24
%mem=175gb
#p wB97XD/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 35.71600 1.25500 29.88400
C  33.26000 2.26900 32.06200
C  37.97500 1.18400 32.63200
C  38.34600 0.78200 27.81300
C  33.67000 1.64600 27.25500
N  35.61600 1.68700 32.10200
C  34.47500 2.10300 32.73400
C  34.68900 2.20300 34.24600
C  36.27600 2.03500 34.33300
C  36.65900 1.61200 32.94500
C  37.02900 3.30200 34.88200
C  33.91200 1.08400 35.07900
C  33.16700 1.39000 36.35600
H  33.84300 2.22100 37.38800
N  37.88700 0.95800 30.15700
C  38.55000 0.78800 31.32200
C  39.89100 0.36100 31.10100
C  40.02700 0.20000 29.68600
C  38.73000 0.72700 29.16600
C  40.80200 0.22400 32.27800
C  41.26400 -0.16900 28.88500
O  41.33400 -0.00500 27.66300
C  42.48400 -0.74400 29.58300
N  36.02000 1.36700 27.77900
C  37.17800 1.00100 27.16900
C  36.98100 0.89100 25.65800
C  35.41500 1.00100 25.51900
C  34.99600 1.35200 26.93300
C  37.78500 1.78800 24.67100
C  34.66500 -0.24600 24.88900
C  34.10400 -1.33400 25.78000
N  33.89900 1.92200 29.64700
C  33.07800 1.91600 28.54800
C  31.77800 2.41500 28.81300
C  31.77700 2.52500 30.20000
C  33.07800 2.19300 30.65800
C  30.56200 2.67700 27.93400
C  31.02900 2.87500 31.32700
O  29.88700 3.29400 31.37200
C  31.95500 2.75800 32.66400
C  32.02000 4.07800 33.27900
O  32.54800 5.09200 32.78100
O  31.38700 4.14700 34.48000
C  31.10000 5.46600 35.05500
H  38.54900 1.11500 33.55800
H  39.09900 0.35600 27.14700
H  32.90400 1.43100 26.50800
H  34.28100 3.15500 34.58800
H  36.55500 1.25900 35.04600
H  37.56800 3.00200 35.78000
H  36.24700 4.05100 35.00800
H  37.83800 3.67300 34.25200
H  34.72200 0.38900 35.30100
H  33.30700 0.45900 34.42300
H  32.77000 0.47000 36.78500
H  32.33500 2.02700 36.05800
H  41.31400 -0.72500 32.12300
H  40.23000 0.11700 33.20000
H  41.52900 1.02400 32.42200
H  42.25000 -1.48700 30.34600
H  43.15800 0.06500 29.86800
H  42.99900 -1.30000 28.79900
H  37.38700 -0.09500 25.43100
H  35.10800 1.75000 24.79000
H  37.12800 2.52600 24.21100
H  38.16400 1.12200 23.89600
H  38.57800 2.24100 25.26500
H  35.37100 -0.71500 24.20500
H  33.81300 0.09600 24.30100
H  34.44300 -2.32400 25.47400
H  33.01400 -1.34700 25.82000
H  34.48200 -1.09900 26.77500
H  30.79100 2.37600 26.91100
H  30.30400 3.73500 27.98700
H  29.72700 2.09100 28.31700
H  31.32600 2.10900 33.27400
H  30.27200 5.47100 35.76400
H  31.01600 6.28300 34.33900
H  31.93600 5.66400 35.72600
Mg -9.05100 41.13500 42.30500
C  -8.30500 37.72200 41.48700
C  -6.39500 41.95900 40.06300
C  -10.13500 44.30600 42.06300
C  -11.94600 40.07300 43.78400
N  -7.59700 39.95400 40.94500
C  -7.46300 38.60400 40.80000
C  -6.30700 38.23500 39.80600
C  -5.51300 39.62600 39.73200
C  -6.60000 40.57800 40.20500
C  -4.30400 39.71300 40.63500
C  -6.79900 37.60800 38.51200
C  -6.03200 36.37000 37.97900
H  -4.96500 36.51100 36.79600
N  -8.29800 42.98100 41.27000
C  -7.13300 43.08000 40.48100
C  -6.94700 44.48300 40.09100
C  -8.03900 45.19400 40.72800
C  -8.91000 44.19100 41.30200
C  -5.89600 45.03300 39.20200
C  -8.19100 46.69400 40.76500
O  -7.26900 47.37100 40.36900
C  -9.41600 47.41600 41.27700
N  -10.80500 41.95700 42.74100
C  -10.93600 43.34200 42.69200
C  -12.32300 43.68900 43.23800
C  -12.78500 42.37500 43.96800
C  -11.75000 41.40700 43.50300
C  -12.48200 45.00600 44.08300
C  -14.25700 41.90100 43.76600
C  -15.09300 41.67300 44.98300
N  -9.92100 39.32400 42.71000
C  -11.08300 39.02100 43.39800
C  -11.16100 37.57500 43.68200
C  -10.08800 37.04400 42.86600
C  -9.41700 38.14200 42.31900
C  -12.14700 36.79800 44.51400
C  -9.42300 35.76800 42.46900
O  -9.72600 34.63000 42.68900
C  -8.33800 36.19100 41.48400
C  -7.03500 35.63000 41.97000
O  -6.57500 36.04000 43.01400
O  -6.38300 34.71100 41.09400
C  -4.99400 34.30500 41.46700
H  -5.55800 42.25300 39.42600
H  -10.61600 45.28500 42.09600
H  -12.77900 39.78900 44.43000
H  -5.57300 37.57400 40.26600
H  -5.28500 39.90000 38.70200
H  -4.64500 39.94800 41.64300
H  -3.63900 40.46300 40.20700
H  -3.69100 38.82300 40.78100
H  -6.68900 38.37600 37.74700
H  -7.87800 37.47800 38.59500
H  -6.83500 35.70800 37.65400
H  -5.48300 35.83500 38.75400
H  -5.12500 45.43600 39.85800
H  -6.29100 45.90400 38.67900
H  -5.42500 44.40000 38.44900
H  -9.41600 48.49600 41.13100
H  -9.43700 47.22500 42.35000
H  -10.40300 47.09200 40.94900
H  -12.86200 43.92500 42.32000
H  -12.52900 42.46800 45.02300
H  -11.51700 45.33200 44.47200
H  -13.22300 44.86700 44.87000
H  -12.90400 45.79900 43.46600
H  -14.25100 40.98200 43.18000
H  -14.72900 42.61100 43.08700
H  -15.22400 40.60100 45.12900
H  -16.08000 42.11100 44.83500
H  -14.72100 42.06100 45.93100
H  -11.73100 35.83800 44.81900
H  -13.03600 36.66900 43.89600
H  -12.46300 37.43800 45.33800
H  -8.66000 35.87700 40.49100
H  -5.04200 33.54700 42.24900
H  -4.44500 35.14400 41.89400
H  -4.44900 33.94400 40.59500


