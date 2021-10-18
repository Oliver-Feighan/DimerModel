%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 16.41100 52.36700 25.56600
C  18.01800 50.88800 28.30000
C  13.92400 53.47100 27.63400
C  14.78500 53.43800 22.89200
C  18.97700 50.98900 23.49400
N  16.01700 52.13900 27.73500
C  16.81400 51.56600 28.70100
C  16.31600 51.93200 30.08700
C  14.87100 52.50500 29.80500
C  15.00600 52.79900 28.27200
C  13.75200 51.51700 30.10900
C  17.38000 52.82700 30.82100
C  16.86600 53.42500 32.16000
H  15.83800 52.61100 32.91800
N  14.67300 53.30500 25.31300
C  13.70900 53.59600 26.27000
C  12.51400 54.21500 25.67000
C  12.80900 54.23100 24.27200
C  14.13200 53.63200 24.11700
C  11.36800 54.71200 26.46900
C  11.91000 54.75300 23.18100
O  12.26900 54.77100 21.97100
C  10.56700 55.35400 23.49600
N  16.86500 52.27500 23.47100
C  15.99700 52.81000 22.61200
C  16.44300 52.50800 21.11400
C  17.81800 51.65600 21.31000
C  17.89300 51.61900 22.83700
C  15.27500 51.72000 20.39100
C  19.11700 52.36100 20.72800
C  19.40800 52.10700 19.22100
N  18.08500 51.15100 25.76500
C  19.11600 50.75400 24.87300
C  20.19800 50.00100 25.64600
C  19.78200 49.96200 26.97300
C  18.52400 50.75100 27.00100
C  21.48500 49.44800 25.11900
C  20.10700 49.55500 28.32000
O  21.06000 48.92400 28.78800
C  19.05800 50.30100 29.28500
C  18.45300 49.37700 30.31700
O  17.48000 48.69900 30.12500
O  19.02600 49.59100 31.59100
C  18.26600 49.12100 32.75900
H  13.19300 53.70900 28.41000
H  14.23600 53.78000 22.01200
H  19.73200 50.45400 22.91400
H  16.24800 51.01000 30.66300
H  14.63800 53.45400 30.28900
H  13.10300 52.02100 30.82500
H  14.20000 50.58000 30.43900
H  13.06200 51.39100 29.27500
H  17.73700 53.44800 30.00000
H  18.21600 52.16100 31.03200
H  16.46200 54.39100 31.85700
H  17.72400 53.56200 32.81800
H  10.55700 54.06200 26.14100
H  11.17100 55.77400 26.32500
H  11.38200 54.60100 27.55400
H  10.71400 56.25300 24.09300
H  10.12800 54.49800 24.01000
H  9.92000 55.56900 22.64600
H  16.62800 53.52500 20.76600
H  17.77900 50.58400 21.11200
H  14.41400 51.50600 21.02400
H  15.64500 50.76000 20.03000
H  14.90700 52.39400 19.61800
H  19.98300 52.23400 21.37700
H  18.92400 53.43100 20.79800
H  19.48300 53.08600 18.74900
H  18.62700 51.45900 18.82300
H  20.36100 51.60000 19.06800
H  22.27100 50.10500 25.49200
H  21.46400 49.41600 24.02900
H  21.78400 48.45300 25.44800
H  19.61000 51.11500 29.75400
H  18.44700 48.04700 32.79300
H  17.22900 49.45300 32.70500
H  18.68800 49.63700 33.62100
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


