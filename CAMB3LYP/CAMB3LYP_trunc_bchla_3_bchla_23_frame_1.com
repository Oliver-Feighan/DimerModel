%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 1.39700 7.67200 26.33400
C  1.72400 9.74100 29.05100
C  2.01300 4.99100 28.38100
C  1.18200 5.57500 23.58000
C  1.22200 10.43500 24.24700
N  1.88100 7.36900 28.50400
C  1.86000 8.39700 29.40300
C  2.08800 7.91400 30.75300
C  2.40000 6.43300 30.60700
C  2.12200 6.20500 29.05700
C  3.81300 5.91300 31.04400
C  0.86900 8.21500 31.69000
C  1.18400 8.50600 33.19800
H  2.44800 7.85400 33.79200
N  1.46600 5.54900 26.07500
C  1.71900 4.62900 27.07800
C  1.70700 3.27200 26.49500
C  1.45400 3.42700 25.06400
C  1.33100 4.87900 24.84700
C  1.96100 1.99300 27.27600
C  1.22600 2.33500 24.04500
O  0.94200 2.65000 22.87600
C  1.34900 0.88300 24.32000
N  1.30400 7.94800 24.20100
C  1.17000 6.94200 23.33200
C  0.89900 7.47700 21.91600
C  0.64800 9.06900 22.12400
C  1.00600 9.14600 23.65000
C  2.06300 7.12000 20.94100
C  -0.80600 9.55500 21.87400
C  -1.00500 10.45600 20.64600
N  1.52400 9.70800 26.54100
C  1.45100 10.73400 25.59600
C  1.56500 12.05000 26.26300
C  1.69000 11.65900 27.61700
C  1.66900 10.24300 27.72900
C  1.49300 13.37100 25.67200
C  1.73800 12.19000 28.99500
O  1.68200 13.31300 29.42100
C  1.80500 10.95700 29.96500
C  2.97400 11.12000 30.84400
O  4.12200 10.98100 30.39200
O  2.68100 11.28600 32.16300
C  3.77100 11.25200 33.13300
H  2.21100 4.15900 29.06000
H  1.15400 4.90100 22.72100
H  1.03400 11.27000 23.56900
H  2.96600 8.42700 31.14700
H  1.62000 5.81700 31.05600
H  4.56300 6.57700 31.47500
H  4.27400 5.43500 30.18000
H  3.77100 5.07400 31.73900
H  0.12500 7.42500 31.58700
H  0.27300 9.01100 31.24400
H  0.33700 8.21300 33.81900
H  1.41300 9.55400 33.38900
H  2.02800 2.15000 28.35200
H  2.91300 1.54300 26.99300
H  1.08700 1.34400 27.22400
H  1.03600 0.37200 23.41000
H  0.74800 0.54700 25.16600
H  2.36800 0.65900 24.63600
H  -0.06900 7.08000 21.61000
H  1.35900 9.68000 21.56800
H  1.72500 6.67300 20.00600
H  2.74700 6.40500 21.39900
H  2.59200 8.06400 20.80800
H  -1.22600 10.10700 22.71400
H  -1.49800 8.72100 21.75500
H  -0.24900 10.22000 19.89700
H  -0.86300 11.48900 20.96200
H  -2.04500 10.34200 20.34300
H  1.26400 14.14800 26.40100
H  0.76200 13.44400 24.86700
H  2.50500 13.65900 25.38700
H  0.85200 11.05500 30.48600
H  4.47400 12.08100 33.04700
H  4.27000 10.28300 33.12600
H  3.36500 11.36600 34.13800
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


