%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 7.14200 57.54600 41.51700
C  6.28900 54.20700 41.28700
C  10.28800 56.72100 40.38900
C  7.81200 60.89700 41.16900
C  3.94100 58.30400 42.53200
N  8.10000 55.73200 40.72500
C  7.55000 54.46800 40.80800
C  8.51900 53.39300 40.32300
C  9.88500 54.17400 40.37000
C  9.39700 55.61300 40.45200
C  10.94100 53.73800 41.48400
C  8.18100 52.71900 38.95500
C  8.31100 51.21100 38.91500
H  9.48100 50.75900 38.05500
N  8.80400 58.60500 40.85800
C  9.98200 58.10200 40.47200
C  10.93200 59.16600 40.19700
C  10.33500 60.40200 40.38100
C  8.90800 60.02900 40.77200
C  12.41400 58.88100 39.65800
C  11.01900 61.80700 40.27800
O  12.22600 61.86700 39.90900
C  10.28000 63.06100 40.55000
N  5.92900 59.42300 41.58100
C  6.46000 60.62400 41.57700
C  5.53200 61.67900 42.14200
C  4.20400 60.94400 42.52900
C  4.70000 59.43900 42.28400
C  6.08200 62.75400 43.13300
C  2.88100 61.24800 41.76800
C  1.64800 60.78800 42.51800
N  5.41300 56.46200 41.84200
C  4.16700 56.89800 42.34900
C  3.34900 55.73600 42.56600
C  4.13900 54.62000 42.19700
C  5.33900 55.15600 41.73800
C  2.01600 55.69000 43.29800
C  4.11700 53.24300 42.14700
O  3.28800 52.39200 42.32200
C  5.50100 52.88200 41.39400
C  5.07100 52.59100 39.95700
O  4.83200 53.42500 39.10500
O  4.99200 51.21300 39.83600
C  4.48800 50.71300 38.56100
H  11.34700 56.50800 40.23000
H  8.09700 61.94400 41.04900
H  2.94700 58.37800 42.97800
H  8.57300 52.63000 41.09900
H  10.29300 54.14500 39.35900
H  11.53300 54.64700 41.59000
H  11.67500 53.00700 41.14400
H  10.53800 53.39000 42.43500
H  8.89800 53.17500 38.27200
H  7.26400 53.03800 38.46000
H  7.37100 50.81600 38.52900
H  8.48000 50.78000 39.90200
H  12.65000 57.83100 39.48000
H  13.15600 59.27900 40.35100
H  12.59000 59.28000 38.66000
H  9.90800 63.04900 41.57400
H  9.55300 63.27700 39.76700
H  11.06400 63.81300 40.63000
H  5.09100 62.35700 41.41100
H  4.06100 60.88100 43.60800
H  7.16900 62.67900 43.08000
H  5.65900 62.68700 44.13600
H  5.85900 63.78200 42.84700
H  2.99200 60.79700 40.78200
H  2.91000 62.32500 41.60300
H  0.83700 61.47100 42.26400
H  1.73700 60.78400 43.60500
H  1.32700 59.82700 42.11600
H  1.67100 54.68100 43.52400
H  1.28900 56.22800 42.69100
H  2.06800 56.22200 44.24800
H  6.04900 52.10100 41.92200
H  3.44600 51.01000 38.44100
H  4.61600 49.63500 38.46900
H  5.02300 51.19800 37.74500
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


