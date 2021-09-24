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
Mg -2.59000 34.01700 26.79400
C  -3.62000 32.41400 29.77300
C  -1.15300 36.45800 28.71900
C  -2.20600 35.82500 24.04600
C  -4.51700 31.70900 25.00200
N  -2.55200 34.47200 29.00300
C  -2.95200 33.59300 30.01300
C  -2.37700 34.05700 31.37400
C  -1.60200 35.37500 31.00700
C  -1.71600 35.46600 29.46000
C  -1.91400 36.61400 31.82800
C  -1.50900 33.04200 32.13400
C  -1.98900 32.76400 33.58500
H  -0.95800 32.99700 34.66700
N  -1.70200 35.86000 26.44900
C  -1.19200 36.74900 27.34900
C  -0.55100 37.80700 26.60500
C  -0.98500 37.71000 25.21400
C  -1.68000 36.43200 25.23200
C  0.22600 39.02700 27.32500
C  -1.03600 38.70200 24.13300
O  -1.40600 38.53400 23.00600
C  -0.57600 40.11000 24.49900
N  -3.46600 33.87800 24.82300
C  -3.08300 34.73100 23.85300
C  -3.52500 34.22200 22.47900
C  -3.92700 32.76700 22.76900
C  -3.94500 32.73700 24.28100
C  -4.74000 35.08300 21.90100
C  -2.77000 31.79500 22.20500
C  -1.44300 31.78800 23.01200
N  -3.75200 32.29300 27.19100
C  -4.44100 31.48200 26.37200
C  -5.05500 30.45500 27.17000
C  -4.72100 30.75200 28.47800
C  -3.96900 31.93000 28.47700
C  -5.92100 29.29000 26.73000
C  -4.97000 30.27500 29.84100
O  -5.63000 29.44100 30.39400
C  -4.20300 31.28400 30.73600
C  -5.22200 31.86200 31.64700
O  -6.22400 32.38000 31.22200
O  -4.94100 31.66200 32.98900
C  -6.05600 31.94700 33.95300
H  -0.52600 37.10900 29.33200
H  -2.00600 36.32000 23.09300
H  -4.95700 30.83400 24.51900
H  -3.20300 34.37600 32.01000
H  -0.55500 35.23200 31.27500
H  -2.84600 36.49700 32.37900
H  -1.95200 37.53200 31.24100
H  -1.13500 36.66900 32.58800
H  -0.45800 33.28400 32.29200
H  -1.47400 32.09100 31.60300
H  -2.46500 31.78700 33.66800
H  -2.77300 33.48500 33.81300
H  1.18600 39.15100 26.82400
H  0.47900 38.88500 28.37500
H  -0.26700 39.99700 27.38900
H  0.46900 39.96700 24.77500
H  -1.19700 40.39100 25.34900
H  -0.59300 40.82400 23.67500
H  -2.71600 34.17900 21.74900
H  -4.90200 32.39200 22.45700
H  -4.66500 35.04300 20.81400
H  -4.79700 36.14200 22.15500
H  -5.62600 34.53900 22.22800
H  -2.53100 31.98600 21.15900
H  -3.07000 30.75600 22.34200
H  -0.61000 31.26900 22.53900
H  -1.63900 31.28300 23.95800
H  -1.17900 32.82600 23.21700
H  -5.41600 28.57400 26.08200
H  -6.72000 29.76700 26.16300
H  -6.40800 28.78800 27.56600
H  -3.52700 30.64600 31.30500
H  -5.48500 32.48600 34.70900
H  -6.48400 31.04500 34.39000
H  -6.92400 32.52400 33.63500


