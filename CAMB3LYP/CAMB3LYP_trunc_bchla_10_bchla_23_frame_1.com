%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 40.74100 8.40300 29.69900
C  42.66000 10.02500 32.01600
C  38.73200 7.21500 32.07700
C  39.17900 6.75200 27.36300
C  43.28300 9.22800 27.34100
N  40.79300 8.45600 31.86800
C  41.50500 9.39200 32.54400
C  40.94800 9.63500 33.92000
C  39.85300 8.49700 34.00200
C  39.82700 7.94800 32.59200
C  40.10000 7.44100 35.10600
C  40.44100 11.18100 34.18500
C  40.87000 11.88200 35.48300
H  41.23900 10.96600 36.69200
N  39.13800 7.34400 29.69200
C  38.39900 6.91700 30.74900
C  37.34700 6.00300 30.37600
C  37.33700 6.02300 29.01400
C  38.65500 6.62800 28.64700
C  36.61000 5.29800 31.50000
C  36.33400 5.36600 28.04500
O  36.44200 5.38700 26.84100
C  35.20700 4.65100 28.73000
N  41.30000 7.91900 27.67200
C  40.36700 7.31400 26.91500
C  40.80200 7.26400 25.44700
C  41.95100 8.32700 25.38700
C  42.23200 8.49800 26.89400
C  41.21600 5.86100 24.82900
C  41.58000 9.72400 24.76100
C  40.48400 10.59900 25.31000
N  42.57200 9.48100 29.57000
C  43.55300 9.66900 28.59100
C  44.67900 10.46300 29.09900
C  44.36400 10.63700 30.45500
C  43.08700 10.02200 30.64800
C  45.85600 11.01300 28.33500
C  44.78800 11.17000 31.66800
O  45.80500 11.75000 32.01500
C  43.61600 11.01200 32.67100
C  44.29700 10.49000 33.87800
O  44.63700 9.31800 34.05700
O  44.61300 11.49900 34.74100
C  45.54500 11.24000 35.84300
H  37.99000 7.00700 32.85100
H  38.68100 6.27000 26.52000
H  43.95700 9.44200 26.50800
H  41.68000 9.37000 34.68300
H  38.87700 8.95900 34.14900
H  40.48500 6.53200 34.64400
H  39.14600 7.20600 35.57800
H  40.83000 7.87400 35.79000
H  39.35300 11.12200 34.15300
H  40.80700 11.84200 33.40000
H  40.04300 12.46600 35.88600
H  41.65400 12.60800 35.26800
H  37.22900 4.66200 32.13300
H  35.89000 4.51000 31.28000
H  36.12800 6.14400 31.99000
H  34.67500 5.32500 29.40100
H  35.72300 3.80100 29.17600
H  34.37800 4.27900 28.12900
H  39.88500 7.56300 24.93800
H  42.83200 7.86300 24.94300
H  41.11200 5.08600 25.58800
H  42.28100 5.83900 24.60100
H  40.70000 5.56000 23.91600
H  41.36300 9.40900 23.74000
H  42.43100 10.40400 24.73900
H  39.92000 10.04700 26.06300
H  39.77700 10.86400 24.52400
H  40.79900 11.50000 25.83600
H  46.78200 10.48700 28.56300
H  45.97200 12.07500 28.55200
H  45.60900 10.77100 27.30100
H  43.07400 11.94700 32.81400
H  46.05900 12.14500 36.16900
H  46.18500 10.38200 35.64300
H  45.04400 10.86200 36.73400
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


