%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 24.45900 -7.80400 45.73000
C  26.58600 -5.04500 44.57700
C  21.88700 -6.07400 44.19000
C  22.84100 -10.55200 45.75100
C  27.52200 -9.37000 46.76800
N  24.34300 -5.76600 44.35600
C  25.30500 -4.83000 44.04500
C  24.68800 -3.74900 43.08400
C  23.16700 -3.94300 43.21100
C  23.11900 -5.33100 43.95600
C  22.33200 -2.87300 44.03000
C  25.23700 -3.82200 41.61800
C  25.48800 -2.49900 40.92700
H  24.41700 -2.03400 39.89100
N  22.60400 -8.23500 45.19100
C  21.64700 -7.40100 44.68000
C  20.38100 -8.15500 44.46900
C  20.67900 -9.47200 44.90600
C  22.09500 -9.47000 45.26500
C  19.08200 -7.63000 43.83600
C  19.81800 -10.73300 44.78900
O  18.67600 -10.59500 44.35900
C  20.22700 -12.02800 45.16700
N  25.20700 -9.77600 45.95500
C  24.19400 -10.69000 46.16700
C  24.84300 -12.02000 46.71400
C  26.38100 -11.70100 46.72900
C  26.35300 -10.16200 46.56900
C  24.30300 -12.52000 48.08200
C  27.32200 -12.42400 45.68300
C  27.76100 -13.86800 45.98200
N  26.63200 -7.37200 45.76700
C  27.69000 -8.05400 46.33400
C  28.89800 -7.26500 46.33900
C  28.50100 -6.05100 45.66300
C  27.12300 -6.13300 45.25800
C  30.21600 -7.74700 46.86200
C  28.99300 -4.84100 45.17500
O  30.13700 -4.41700 45.08200
C  27.84800 -4.18900 44.34700
C  27.67100 -2.81900 44.95900
O  27.62100 -1.78000 44.29600
O  27.39300 -2.80300 46.38500
C  26.88100 -1.51100 46.83900
H  21.04400 -5.51700 43.77700
H  22.37000 -11.50900 45.98500
H  28.34500 -9.94300 47.20000
H  24.83100 -2.76600 43.53500
H  22.69900 -4.12200 42.24300
H  21.79800 -2.30300 43.27000
H  22.98800 -2.17400 44.54800
H  21.67600 -3.27700 44.80100
H  24.55300 -4.27200 40.89800
H  26.19500 -4.34200 41.59800
H  26.49400 -2.54600 40.51000
H  25.50300 -1.72600 41.69600
H  18.22200 -8.02700 44.37600
H  19.15100 -7.95200 42.79700
H  19.01200 -6.54700 43.92800
H  20.35400 -12.03800 46.25000
H  21.09800 -12.40800 44.63300
H  19.40000 -12.62300 44.77800
H  24.71500 -12.87500 46.05100
H  26.80800 -11.85300 47.72000
H  25.17900 -12.58000 48.72800
H  23.87600 -13.51200 47.93800
H  23.57100 -11.76000 48.35800
H  28.20300 -11.78200 45.65400
H  26.82800 -12.52100 44.71600
H  27.32600 -14.00800 46.97200
H  28.83900 -13.92900 46.13800
H  27.42200 -14.64300 45.29500
H  30.43000 -8.65700 46.30300
H  30.17700 -8.08000 47.89900
H  31.04400 -7.08900 46.59800
H  28.22500 -4.06500 43.33200
H  26.34600 -1.90700 47.70300
H  26.35400 -0.97600 46.04900
H  27.73900 -0.96400 47.23000
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


