%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 44.80500 2.87200 46.67800
C  42.56400 5.57100 46.81800
C  42.13800 0.79200 46.35700
C  47.01300 0.52400 46.12300
C  47.34200 5.37800 46.17800
N  42.64100 3.17100 46.66900
C  41.94200 4.34400 46.71700
C  40.48300 4.04900 46.33800
C  40.34000 2.59000 46.75700
C  41.78900 2.14300 46.48400
C  39.96100 2.42600 48.24700
C  40.23600 4.25700 44.79800
C  39.26900 5.45400 44.33400
H  39.83900 6.57400 43.44700
N  44.60700 0.93700 46.23600
C  43.41900 0.24600 46.16900
C  43.72000 -1.14400 46.02500
C  45.14200 -1.32500 46.06900
C  45.67600 0.06600 46.14000
C  42.62100 -2.25300 45.89000
C  45.87200 -2.71700 46.06400
O  45.23700 -3.72500 46.11700
C  47.31400 -2.74100 46.27700
N  46.83800 2.95400 45.90000
C  47.54500 1.82700 46.15300
C  49.01900 2.09000 46.41200
C  49.17500 3.56200 46.06200
C  47.67600 4.02200 46.03300
C  49.44600 1.78300 47.89500
C  49.96600 3.97900 44.73400
C  51.33300 4.60700 44.78400
N  45.00300 4.99200 46.72600
C  46.07700 5.86500 46.51000
C  45.65700 7.24900 46.55200
C  44.29400 7.16900 46.69900
C  43.89900 5.79300 46.77900
C  46.53700 8.41600 46.37700
C  43.08200 7.96600 46.81100
O  42.95000 9.19400 46.88400
C  41.94500 6.98400 47.00000
C  41.32800 7.26400 48.36600
O  42.00600 7.20100 49.38900
O  40.00100 7.70700 48.25300
C  39.36700 8.02300 49.50500
H  41.26200 0.14300 46.41400
H  47.71200 -0.30200 46.26300
H  48.13300 6.13000 46.13700
H  39.87300 4.80700 46.83000
H  39.59400 2.12600 46.11300
H  39.12500 1.73800 48.12600
H  39.66800 3.42600 48.56800
H  40.73100 2.03300 48.91100
H  39.76200 3.35000 44.42300
H  41.21000 4.30800 44.31200
H  38.98200 5.96200 45.25500
H  38.33700 5.13200 43.87100
H  42.71600 -2.97100 45.07600
H  41.68000 -1.71700 45.76200
H  42.51800 -2.85500 46.79300
H  47.88700 -2.17400 45.54300
H  47.64600 -3.77900 46.25800
H  47.54400 -2.49600 47.31300
H  49.67000 1.44400 45.82200
H  49.62700 4.10200 46.89400
H  49.91800 0.80800 48.01500
H  48.60400 1.96800 48.56200
H  50.20200 2.52100 48.16700
H  49.41200 4.48800 43.94600
H  50.09500 2.95900 44.37100
H  52.11800 3.97300 44.37300
H  51.60000 4.81700 45.82000
H  51.22600 5.52900 44.21300
H  47.39900 8.09600 45.79100
H  46.88900 8.66700 47.37800
H  45.99500 9.26600 45.96500
H  41.25100 7.27600 46.21300
H  39.04600 7.15600 50.08300
H  38.48400 8.64500 49.36100
H  39.98200 8.61300 50.18400
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


