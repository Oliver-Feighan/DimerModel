%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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


