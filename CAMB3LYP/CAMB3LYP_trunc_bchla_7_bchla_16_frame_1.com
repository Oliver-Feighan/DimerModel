%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 26.07000 0.32300 29.17500
C  27.93000 0.00300 32.17200
C  23.23700 0.37300 31.22200
C  24.20000 0.73500 26.50800
C  28.89900 -0.55900 27.35600
N  25.68800 0.23300 31.42600
C  26.58700 0.25300 32.44900
C  25.99300 0.25100 33.86800
C  24.43500 0.24800 33.45400
C  24.48700 0.27000 31.95500
C  23.63400 -1.02300 33.94300
C  26.39900 1.50900 34.70300
C  27.30600 1.47800 35.92700
H  27.07500 2.41400 37.10600
N  24.00400 0.60600 28.87400
C  22.98600 0.52400 29.83500
C  21.72400 0.61300 29.13500
C  21.97300 0.88000 27.79100
C  23.44300 0.80400 27.68300
C  20.46800 0.64800 29.93700
C  21.02800 1.08200 26.60700
O  21.39200 1.28100 25.44400
C  19.58500 1.18400 26.91000
N  26.42900 0.02500 27.20000
C  25.57000 0.46200 26.23000
C  26.30100 0.60400 24.85200
C  27.70600 0.10100 25.15700
C  27.68000 -0.24400 26.67200
C  25.55900 -0.10800 23.67400
C  28.73300 1.18300 24.79000
C  29.96300 0.77300 24.06500
N  28.00500 -0.28200 29.61000
C  29.05600 -0.58000 28.78000
C  30.24500 -0.86800 29.52500
C  29.92400 -0.56000 30.93700
C  28.52000 -0.19800 30.89800
C  31.58600 -1.33700 28.93500
C  30.40700 -0.48700 32.33900
O  31.55800 -0.58600 32.79700
C  29.11500 -0.11100 33.17500
C  29.05800 -1.28400 34.10300
O  28.53100 -2.37500 33.92200
O  29.56000 -0.97800 35.26800
C  29.75400 -2.03800 36.28700
H  22.46700 0.52900 31.98000
H  23.69600 0.81600 25.54300
H  29.79300 -0.76400 26.76400
H  26.26300 -0.69600 34.33700
H  23.82700 1.11200 33.72400
H  24.21800 -1.88400 34.26700
H  22.99700 -1.39000 33.13900
H  22.99100 -0.72400 34.77100
H  25.46900 1.94600 35.06700
H  26.90100 2.23700 34.06700
H  28.23400 1.84700 35.49000
H  27.27900 0.42100 36.19200
H  19.92700 1.58400 29.80500
H  20.67800 0.53600 31.00100
H  19.76300 -0.13800 29.66800
H  19.43300 1.90200 27.71600
H  19.18900 0.27500 27.36300
H  18.97400 1.53500 26.07900
H  26.25300 1.68100 24.69200
H  27.96600 -0.74800 24.52500
H  24.56100 -0.46700 23.92500
H  26.08600 -0.97800 23.28400
H  25.56500 0.53500 22.79400
H  29.00400 1.58600 25.76500
H  28.34900 1.97700 24.15000
H  30.78000 1.14300 24.68400
H  29.93300 1.13400 23.03700
H  30.14200 -0.30200 24.03500
H  31.60000 -2.42600 28.91000
H  32.38900 -0.93100 29.55000
H  31.71700 -1.01400 27.90200
H  29.35100 0.80400 33.71800
H  28.83300 -2.53300 36.59300
H  30.27100 -1.70300 37.18600
H  30.41000 -2.86500 36.01700
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


