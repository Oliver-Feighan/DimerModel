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


