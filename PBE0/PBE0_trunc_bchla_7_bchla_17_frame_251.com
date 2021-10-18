%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 25.94600 0.51000 29.41700
C  27.75500 0.04800 32.42200
C  23.08400 0.98300 31.28700
C  24.18800 0.65700 26.52400
C  28.83500 -0.18800 27.55800
N  25.56000 0.39300 31.59400
C  26.46800 0.41400 32.65000
C  25.73500 0.61600 34.03900
C  24.24300 0.62700 33.61200
C  24.25900 0.61600 32.06700
C  23.30400 -0.50900 34.16600
C  26.25400 1.80200 34.95100
C  26.76400 1.39400 36.30000
H  26.45300 2.38700 37.48900
N  23.95800 0.76900 29.03000
C  22.92100 1.01800 29.90600
C  21.67800 1.14300 29.16100
C  21.94400 1.05800 27.76500
C  23.43600 0.83800 27.70800
C  20.32000 1.49500 29.80100
C  21.00100 1.19500 26.57400
O  21.36800 1.18900 25.41700
C  19.44600 1.46900 26.69600
N  26.43900 0.16900 27.31600
C  25.54800 0.29300 26.34100
C  26.18900 0.06900 24.92500
C  27.73100 -0.01000 25.26000
C  27.69500 0.01500 26.81800
C  25.54500 -1.13400 24.17700
C  28.51900 1.10000 24.59300
C  29.85400 0.64900 24.00300
N  27.89800 -0.11200 29.79400
C  28.99000 -0.26300 28.96100
C  30.08800 -0.69100 29.77200
C  29.69800 -0.55100 31.12700
C  28.32100 -0.19300 31.11800
C  31.42300 -1.10700 29.28200
C  30.17600 -0.59300 32.52000
O  31.23000 -0.96300 33.02900
C  28.81000 -0.12800 33.44600
C  28.50200 -1.15100 34.39600
O  27.57900 -1.93000 34.29600
O  29.31100 -1.04800 35.50500
C  29.06500 -2.12600 36.49500
H  22.18200 1.26000 31.83700
H  23.66500 0.64600 25.56600
H  29.70300 -0.30600 26.90500
H  25.82900 -0.35200 34.53000
H  23.76600 1.56200 33.90600
H  22.35300 -0.03400 34.40900
H  23.77900 -0.95200 35.04100
H  23.10600 -1.31400 33.45800
H  25.42400 2.50500 35.02300
H  27.09400 2.25600 34.42500
H  27.85000 1.29700 36.26800
H  26.30600 0.41800 36.45800
H  20.08400 2.49100 29.42800
H  20.46000 1.63600 30.87300
H  19.49000 0.83300 29.55100
H  19.31200 2.38400 27.27400
H  18.87900 0.57700 26.96400
H  19.10200 1.53100 25.66400
H  26.08400 0.92600 24.26000
H  28.04800 -1.00900 24.96000
H  26.34500 -1.78700 23.82700
H  24.86300 -0.83400 23.38100
H  24.93000 -1.70700 24.86900
H  28.70300 1.93600 25.26700
H  27.98200 1.61200 23.79300
H  30.10400 -0.36600 24.31200
H  30.60900 1.32100 24.41200
H  29.90400 0.79600 22.92400
H  31.58400 -1.10700 28.20400
H  31.68500 -2.00300 29.84600
H  32.04200 -0.24000 29.51400
H  29.03700 0.80200 33.96800
H  29.24700 -3.09800 36.03600
H  28.05200 -1.96800 36.86500
H  29.74500 -2.08500 37.34500
Mg 29.10200 58.37000 41.70600
C  26.21100 56.85000 40.52500
C  31.00300 55.92400 39.98500
C  31.78200 60.36200 41.82300
C  27.09800 61.07100 42.76300
N  28.62400 56.74200 40.07100
C  27.36300 56.19900 40.01300
C  27.34200 54.94600 39.19500
C  28.87800 54.50900 39.24100
C  29.59900 55.79900 39.77000
C  29.21600 53.21300 40.13700
C  26.66600 55.21800 37.75400
C  26.06800 54.05700 37.13700
H  26.60500 53.76400 35.72000
N  31.13700 58.12700 41.18100
C  31.72600 57.03700 40.48100
C  33.15300 57.29600 40.39000
C  33.42500 58.57400 40.93600
C  32.09000 59.05300 41.42900
C  34.07200 56.38900 39.66000
C  34.85200 59.26600 40.92400
O  35.77500 58.62800 40.42000
C  35.15900 60.50100 41.76400
N  29.32700 60.43000 42.12600
C  30.57900 60.94600 42.21600
C  30.58100 62.56100 42.56700
C  29.07200 62.72600 43.08400
C  28.42300 61.36300 42.68500
C  31.60200 63.02200 43.67500
C  28.34100 63.90500 42.43900
C  27.75100 65.04300 43.30600
N  27.08200 58.74100 41.92000
C  26.43600 59.84700 42.39100
C  25.02600 59.57700 42.45300
C  24.91200 58.38300 41.64600
C  26.21500 57.96500 41.32600
C  23.95900 60.50800 42.95200
C  23.95200 57.50800 41.00600
O  22.77000 57.42900 41.15100
C  24.75000 56.37400 40.25200
C  24.55300 54.95400 40.82700
O  24.66200 54.66400 41.99900
O  24.13200 54.08300 39.80100
C  23.73600 52.76500 40.19600
H  31.53000 55.10300 39.49600
H  32.55200 61.12000 41.66600
H  26.51800 61.92400 43.12000
H  26.89400 54.10800 39.73000
H  29.20300 54.31600 38.21800
H  29.65600 52.49600 39.44400
H  28.33800 52.79800 40.63100
H  29.96300 53.53300 40.86400
H  27.50600 55.54800 37.14300
H  26.02700 56.07100 37.98400
H  24.99800 54.25700 37.09500
H  26.28800 53.16600 37.72500
H  33.60300 55.42900 39.44500
H  34.92200 56.13300 40.29300
H  34.37800 56.81900 38.70600
H  34.62000 61.39400 41.44600
H  36.23600 60.66600 41.71900
H  34.91600 60.25600 42.79800
H  30.70200 63.02600 41.58900
H  28.95600 62.73400 44.16800
H  32.36800 63.65500 43.22800
H  32.08100 62.14600 44.11400
H  31.02600 63.57600 44.41700
H  27.57000 63.56900 41.74600
H  28.92000 64.44500 41.69000
H  26.66500 65.11700 43.36900
H  28.07200 65.96200 42.81500
H  28.11700 64.98900 44.33200
H  23.07500 59.90100 43.14900
H  23.49000 61.12600 42.18700
H  24.26100 61.15900 43.77200
H  24.32100 56.40300 39.25000
H  23.95500 52.21700 39.27900
H  22.65500 52.67800 40.30400
H  24.24700 52.46200 41.10900


