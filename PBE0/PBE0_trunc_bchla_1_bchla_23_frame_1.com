%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg -1.96800 16.99100 27.07800
C  -2.18500 14.77300 29.87900
C  -3.22300 19.43300 29.10400
C  -2.21300 18.84100 24.51400
C  -1.82200 14.16700 25.04100
N  -2.58500 16.98200 29.20400
C  -2.56200 16.08300 30.18300
C  -2.90400 16.68200 31.58000
C  -3.74100 17.95500 31.11000
C  -3.12300 18.15900 29.73800
C  -5.23100 17.67600 30.94100
C  -1.62300 17.10100 32.37700
C  -1.84400 17.38800 33.86500
H  -0.80000 17.02400 34.81100
N  -2.15300 18.95200 26.92600
C  -2.67100 19.79100 27.85500
C  -2.58800 21.14700 27.27100
C  -2.15500 21.00600 25.85300
C  -2.18200 19.57400 25.69800
C  -2.85500 22.39700 28.05200
C  -1.78700 22.00400 24.80100
O  -1.66500 21.68100 23.57000
C  -1.65800 23.40500 25.24100
N  -2.23300 16.57800 25.04900
C  -2.24700 17.56000 24.14100
C  -2.16400 17.01500 22.74600
C  -1.84400 15.45800 22.89800
C  -1.88200 15.35900 24.44300
C  -3.41900 17.19600 21.88100
C  -0.59500 14.92100 22.15500
C  0.65800 15.50100 22.77200
N  -1.81800 14.84700 27.38100
C  -1.79500 13.85500 26.39000
C  -1.79500 12.55300 27.08200
C  -1.83500 12.88000 28.43600
C  -1.88500 14.29600 28.58000
C  -1.81900 11.27600 26.45500
C  -1.87200 12.36200 29.79400
O  -1.76000 11.24400 30.29900
C  -2.22300 13.56300 30.77900
C  -1.25100 13.50900 31.92600
O  -0.03700 13.73600 31.74700
O  -1.93000 13.43500 33.08100
C  -1.02500 13.40800 34.26100
H  -3.82700 20.21100 29.57500
H  -2.18600 19.44400 23.60400
H  -1.69200 13.24400 24.47300
H  -3.52200 16.04200 32.21000
H  -3.55300 18.82200 31.74300
H  -5.45200 16.61400 31.04700
H  -5.47900 18.01800 29.93600
H  -5.82000 18.29200 31.62100
H  -1.29900 18.08700 32.04400
H  -0.83800 16.35200 32.26700
H  -2.81200 17.01700 34.20300
H  -1.90600 18.46800 33.99900
H  -3.31200 22.17700 29.01600
H  -3.61900 22.95400 27.51000
H  -1.96100 23.01100 28.16200
H  -2.60300 23.73500 25.67300
H  -1.30300 23.96100 24.37400
H  -0.92000 23.45700 26.04200
H  -1.36000 17.49200 22.18700
H  -2.66000 14.83500 22.53100
H  -3.23500 17.73100 20.94900
H  -4.17200 17.70500 22.48200
H  -3.93500 16.27100 21.62500
H  -0.61200 15.23200 21.11100
H  -0.51600 13.83400 22.17600
H  0.74600 16.46900 22.28000
H  1.51900 14.87000 22.54900
H  0.52700 15.68800 23.83800
H  -1.10100 10.65800 26.99400
H  -1.56400 11.34900 25.39800
H  -2.85900 10.99200 26.61500
H  -3.25100 13.41000 31.10600
H  -1.27000 12.51300 34.83300
H  -1.13500 14.37200 34.75800
H  0.03200 13.23700 34.05400
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


