%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 46.26900 25.13700 28.87000
C  47.26400 27.26100 31.45200
C  45.06800 22.98000 31.10200
C  46.59300 22.69600 26.49900
C  47.58600 27.40300 26.57500
N  46.02800 25.20400 31.05400
C  46.48600 26.15500 31.85900
C  45.98900 25.86800 33.32200
C  45.39300 24.44600 33.19000
C  45.51200 24.15500 31.69700
C  46.24000 23.30800 33.94900
C  44.92100 26.90700 33.89400
C  45.18000 27.25800 35.31900
H  44.05100 27.15600 36.33300
N  45.85600 23.12600 28.82000
C  45.33200 22.40200 29.85800
C  44.87200 21.10600 29.34000
C  45.32500 21.03600 27.99300
C  46.00600 22.31800 27.74700
C  44.15900 20.03600 30.22700
C  45.21800 19.84800 26.99000
O  45.63600 19.88000 25.83800
C  44.54500 18.61000 27.44900
N  46.98900 25.03900 26.90500
C  47.10000 23.90200 26.09900
C  47.90400 24.25000 24.81900
C  47.61600 25.80400 24.63200
C  47.36900 26.07600 26.14300
C  49.38900 23.83600 24.90100
C  46.42100 26.18700 23.80600
C  45.03400 26.12100 24.37300
N  47.05200 27.04800 28.89500
C  47.58500 27.84200 27.90600
C  48.07800 29.11700 28.43200
C  48.07000 28.90600 29.88400
C  47.43900 27.61600 30.07900
C  48.67900 30.18600 27.67300
C  48.43700 29.46600 31.17900
O  49.02500 30.43400 31.48700
C  47.84100 28.42700 32.26500
C  48.81200 28.03800 33.35100
O  49.55700 27.06300 33.31900
O  48.64800 28.91100 34.42600
C  49.48500 28.68900 35.61900
H  44.47200 22.24900 31.65300
H  46.68200 21.86300 25.79800
H  47.89200 28.10500 25.79700
H  46.81500 25.79600 34.03000
H  44.37900 24.31000 33.56600
H  45.81000 23.21000 34.94600
H  47.29000 23.54400 34.12000
H  46.10500 22.37300 33.40400
H  43.90600 26.51300 33.84000
H  44.83700 27.84500 33.34500
H  45.39600 28.32500 35.27200
H  46.03100 26.79100 35.81400
H  44.95600 19.41500 30.63500
H  43.40600 19.48400 29.66400
H  43.77200 20.42600 31.16900
H  44.70200 18.21100 28.45100
H  44.79100 17.78900 26.77600
H  43.48800 18.66000 27.18600
H  47.38400 23.69600 24.03800
H  48.47400 26.33800 24.22400
H  50.12500 24.61500 24.70200
H  49.61300 22.90600 24.37800
H  49.63200 23.53200 25.91900
H  46.32100 25.62700 22.87700
H  46.56100 27.21900 23.48400
H  44.87400 27.14600 24.70700
H  44.90600 25.48200 25.24700
H  44.33600 25.84200 23.58300
H  48.36400 31.10800 28.16200
H  48.25000 30.14800 26.67200
H  49.76500 30.10200 27.69900
H  47.00200 28.92500 32.75100
H  50.51900 28.64400 35.27600
H  49.25400 27.73000 36.08200
H  49.44300 29.39800 36.44600
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


