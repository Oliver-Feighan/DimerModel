%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 52.53200 24.27200 45.27800
C  49.69500 26.14500 43.97100
C  51.04500 21.43900 44.03200
C  55.55800 22.86200 45.08400
C  53.78600 27.15600 46.56500
N  50.65200 23.90000 43.94600
C  49.60000 24.82500 43.60600
C  48.41500 24.05300 43.01900
C  48.70000 22.56700 43.38700
C  50.22600 22.57000 43.84000
C  47.76900 22.00700 44.56400
C  48.04400 24.34300 41.59400
C  49.09500 25.06000 40.71400
H  48.96900 25.39700 39.21700
N  53.20200 22.33000 44.76700
C  52.42300 21.31600 44.25500
C  53.27600 20.09300 44.13100
C  54.57100 20.56600 44.32700
C  54.53600 21.98000 44.69300
C  52.68600 18.78600 43.73500
C  55.85700 19.70500 44.16100
O  55.72900 18.56500 43.75100
C  57.24900 20.21900 44.32400
N  54.40600 25.03400 45.62800
C  55.48900 24.15600 45.60300
C  56.74800 24.82300 46.15700
C  56.16900 26.07400 46.91200
C  54.72900 26.07700 46.39400
C  57.69100 23.98800 47.06100
C  56.92500 27.42300 46.68900
C  57.79700 27.94700 47.92000
N  51.87600 26.25900 45.36900
C  52.48600 27.32200 46.01900
C  51.58600 28.44900 46.02500
C  50.55800 28.02100 45.17400
C  50.79100 26.66600 44.83200
C  51.91200 29.79900 46.68700
C  49.30500 28.37700 44.58200
O  48.71800 29.46600 44.74300
C  48.58400 27.14400 43.87700
C  48.23600 27.64400 42.53200
O  48.98500 28.23500 41.74700
O  46.97500 27.21700 42.27700
C  46.32700 27.76200 41.04300
H  50.57900 20.54200 43.62000
H  56.53800 22.38300 45.02800
H  54.16100 28.03300 47.09600
H  47.65100 24.44400 43.69100
H  48.71700 21.90900 42.51800
H  47.15300 21.26400 44.05900
H  47.21300 22.83200 45.01000
H  48.44700 21.51600 45.26300
H  47.08700 24.86300 41.63300
H  47.84600 23.33500 41.22800
H  49.99900 24.46000 40.82300
H  49.29400 26.05200 41.12000
H  53.06600 18.57900 42.73500
H  51.60100 18.88400 43.75400
H  52.95100 18.03100 44.47600
H  57.60500 20.97300 43.62300
H  58.00300 19.45100 44.15400
H  57.23000 20.50700 45.37500
H  57.21200 25.04800 45.19600
H  56.15200 25.77400 47.96000
H  58.56700 23.60300 46.53800
H  57.16500 23.10200 47.41500
H  57.94800 24.60600 47.92200
H  56.28400 28.25100 46.38500
H  57.54500 27.23200 45.81300
H  58.71300 28.26700 47.42400
H  58.01300 27.19900 48.68300
H  57.20100 28.70700 48.42400
H  51.05000 30.36300 47.04400
H  52.35800 30.33200 45.84800
H  52.47600 29.64100 47.60600
H  47.71400 26.85200 44.46700
H  45.32600 27.34000 40.94700
H  46.88700 27.48700 40.14900
H  46.32300 28.83200 41.24600
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


