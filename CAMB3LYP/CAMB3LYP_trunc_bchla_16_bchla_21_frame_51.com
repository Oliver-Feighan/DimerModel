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
Mg 16.41100 52.36700 25.56600
C  18.01800 50.88800 28.30000
C  13.92400 53.47100 27.63400
C  14.78500 53.43800 22.89200
C  18.97700 50.98900 23.49400
N  16.01700 52.13900 27.73500
C  16.81400 51.56600 28.70100
C  16.31600 51.93200 30.08700
C  14.87100 52.50500 29.80500
C  15.00600 52.79900 28.27200
C  13.75200 51.51700 30.10900
C  17.38000 52.82700 30.82100
C  16.86600 53.42500 32.16000
H  15.83800 52.61100 32.91800
N  14.67300 53.30500 25.31300
C  13.70900 53.59600 26.27000
C  12.51400 54.21500 25.67000
C  12.80900 54.23100 24.27200
C  14.13200 53.63200 24.11700
C  11.36800 54.71200 26.46900
C  11.91000 54.75300 23.18100
O  12.26900 54.77100 21.97100
C  10.56700 55.35400 23.49600
N  16.86500 52.27500 23.47100
C  15.99700 52.81000 22.61200
C  16.44300 52.50800 21.11400
C  17.81800 51.65600 21.31000
C  17.89300 51.61900 22.83700
C  15.27500 51.72000 20.39100
C  19.11700 52.36100 20.72800
C  19.40800 52.10700 19.22100
N  18.08500 51.15100 25.76500
C  19.11600 50.75400 24.87300
C  20.19800 50.00100 25.64600
C  19.78200 49.96200 26.97300
C  18.52400 50.75100 27.00100
C  21.48500 49.44800 25.11900
C  20.10700 49.55500 28.32000
O  21.06000 48.92400 28.78800
C  19.05800 50.30100 29.28500
C  18.45300 49.37700 30.31700
O  17.48000 48.69900 30.12500
O  19.02600 49.59100 31.59100
C  18.26600 49.12100 32.75900
H  13.19300 53.70900 28.41000
H  14.23600 53.78000 22.01200
H  19.73200 50.45400 22.91400
H  16.24800 51.01000 30.66300
H  14.63800 53.45400 30.28900
H  13.10300 52.02100 30.82500
H  14.20000 50.58000 30.43900
H  13.06200 51.39100 29.27500
H  17.73700 53.44800 30.00000
H  18.21600 52.16100 31.03200
H  16.46200 54.39100 31.85700
H  17.72400 53.56200 32.81800
H  10.55700 54.06200 26.14100
H  11.17100 55.77400 26.32500
H  11.38200 54.60100 27.55400
H  10.71400 56.25300 24.09300
H  10.12800 54.49800 24.01000
H  9.92000 55.56900 22.64600
H  16.62800 53.52500 20.76600
H  17.77900 50.58400 21.11200
H  14.41400 51.50600 21.02400
H  15.64500 50.76000 20.03000
H  14.90700 52.39400 19.61800
H  19.98300 52.23400 21.37700
H  18.92400 53.43100 20.79800
H  19.48300 53.08600 18.74900
H  18.62700 51.45900 18.82300
H  20.36100 51.60000 19.06800
H  22.27100 50.10500 25.49200
H  21.46400 49.41600 24.02900
H  21.78400 48.45300 25.44800
H  19.61000 51.11500 29.75400
H  18.44700 48.04700 32.79300
H  17.22900 49.45300 32.70500
H  18.68800 49.63700 33.62100


