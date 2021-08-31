%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 41.12800 41.85700 26.98600
C  40.12700 44.17900 29.64000
C  41.31000 39.41000 29.48700
C  42.36200 39.82100 24.75900
C  40.70600 44.36300 24.68300
N  40.96700 41.85000 29.29200
C  40.42100 42.87500 30.17300
C  40.33700 42.44400 31.65300
C  40.94600 40.93000 31.52100
C  41.11600 40.67200 30.01100
C  42.22400 40.71400 32.30800
C  38.87900 42.58800 32.18500
C  38.60800 41.92700 33.50300
H  37.80600 42.65600 34.55800
N  41.79000 39.86100 27.08100
C  41.78200 39.01000 28.19200
C  42.26600 37.68700 27.80200
C  42.58400 37.77600 26.38500
C  42.25600 39.19200 25.97700
C  42.31400 36.43000 28.74600
C  43.16500 36.78400 25.44300
O  43.35700 36.97800 24.25000
C  43.44600 35.41200 25.89800
N  41.39000 42.02100 25.04500
C  42.03300 41.09600 24.33000
C  42.18900 41.49200 22.84300
C  41.28100 42.70200 22.76900
C  41.15000 43.12300 24.28400
C  43.62900 41.80100 22.45000
C  39.92300 42.66800 22.11700
C  38.93300 41.74000 22.78800
N  40.64100 43.87900 27.09600
C  40.41100 44.75400 26.01900
C  40.01000 46.03800 26.53900
C  39.83600 45.76900 27.92100
C  40.24100 44.45600 28.22800
C  39.64600 47.28300 25.78400
C  39.34100 46.43200 29.16700
O  38.79000 47.54000 29.34100
C  39.55800 45.41800 30.28800
C  40.41800 46.03500 31.39300
O  41.61900 46.43600 31.36700
O  39.69400 45.95900 32.54600
C  40.38200 46.43100 33.76900
H  41.17000 38.61400 30.22000
H  42.83200 39.17800 24.01200
H  40.57600 45.21400 24.01100
H  40.93200 43.17900 32.19500
H  40.17200 40.21500 31.80000
H  42.06000 39.90200 33.01700
H  42.40500 41.63800 32.85600
H  43.16100 40.69100 31.75200
H  38.23700 42.28000 31.35900
H  38.65800 43.64900 32.30200
H  39.47300 41.65400 34.10700
H  37.99100 41.04400 33.33900
H  41.63300 35.71400 28.28500
H  41.97100 36.63500 29.76000
H  43.34800 36.08500 28.73700
H  43.92500 35.36100 26.87600
H  44.12500 35.02700 25.13800
H  42.45600 34.98600 26.05800
H  41.84700 40.64300 22.25100
H  41.73700 43.52300 22.21500
H  43.81400 41.38100 21.46200
H  44.33400 41.35700 23.15200
H  43.75500 42.87300 22.29400
H  40.06700 42.45700 21.05700
H  39.41200 43.62400 22.23000
H  38.39800 41.12300 22.06700
H  38.14200 42.27600 23.31300
H  39.38000 41.04100 23.49500
H  38.67000 47.63800 26.11400
H  39.60000 47.28100 24.69500
H  40.39300 48.03500 26.03900
H  38.59300 45.21900 30.75500
H  40.75500 45.57700 34.33400
H  39.56300 47.00000 34.20900
H  41.15300 47.17200 33.55400
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


