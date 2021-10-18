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
Mg -9.26600 18.27700 42.96500
C  -5.81100 17.58000 42.98100
C  -8.47400 21.42900 42.13900
C  -12.47500 18.67000 42.48700
C  -9.75200 14.74100 43.14100
N  -7.34400 19.34600 42.56200
C  -6.06200 18.93000 42.69300
C  -5.04000 20.02200 42.45000
C  -5.98500 21.29000 42.45900
C  -7.37700 20.70000 42.36700
C  -5.82600 22.23200 43.72300
C  -4.08200 20.00700 41.18000
C  -4.36600 19.09100 39.96000
H  -3.96500 19.61900 38.53200
N  -10.31800 19.93500 42.50700
C  -9.83700 21.10100 42.09700
C  -10.97100 21.97200 41.68400
C  -12.19800 21.17200 41.95900
C  -11.67600 19.85500 42.32500
C  -10.82500 23.41900 41.29700
C  -13.61600 21.63300 41.74500
O  -13.81100 22.79100 41.38400
C  -14.89900 20.73300 41.88300
N  -10.94800 16.84400 42.71500
C  -12.12100 17.37600 42.70200
C  -13.18200 16.19300 42.82000
C  -12.25700 14.97500 42.77500
C  -10.85500 15.54600 42.86600
C  -14.16200 16.18100 44.00500
C  -12.30000 14.18000 41.36900
C  -12.88700 12.75700 41.40300
N  -8.10600 16.49700 42.98200
C  -8.41600 15.17700 43.17700
C  -7.19300 14.39200 43.48700
C  -6.09500 15.27800 43.29700
C  -6.75800 16.56600 43.10100
C  -7.09800 12.86900 43.72900
C  -4.68100 15.45100 43.26700
O  -3.71100 14.67400 43.28400
C  -4.42500 16.97400 43.13600
C  -3.61000 17.39300 44.36300
O  -2.45500 17.69200 44.51200
O  -4.50700 17.60300 45.41600
C  -3.95600 18.13500 46.70100
H  -8.18700 22.44600 41.86300
H  -13.54300 18.86100 42.35700
H  -9.93500 13.71100 43.45300
H  -4.44700 20.04900 43.36400
H  -5.75600 21.91400 41.59500
H  -5.11000 21.81600 44.43100
H  -6.79600 22.11100 44.20400
H  -5.58800 23.24500 43.39800
H  -3.08000 19.78600 41.54800
H  -3.99700 21.02300 40.79300
H  -5.44200 18.91700 39.97000
H  -3.77000 18.17900 39.94600
H  -11.77100 23.96000 41.28900
H  -10.25700 23.43000 40.36600
H  -10.11700 23.80900 42.02900
H  -15.08400 20.15000 42.78600
H  -14.85100 19.96100 41.11500
H  -15.81700 21.30800 41.76500
H  -13.87700 16.19700 41.98100
H  -12.39400 14.18200 43.51100
H  -14.02600 17.04700 44.65300
H  -14.09900 15.25100 44.57100
H  -15.20700 16.25900 43.70500
H  -11.30400 14.05300 40.94700
H  -12.88100 14.81200 40.69700
H  -12.44900 12.05100 40.69700
H  -13.92900 12.87300 41.10500
H  -12.71100 12.44000 42.43100
H  -7.98800 12.30200 43.45800
H  -6.94100 12.78400 44.80400
H  -6.18900 12.50600 43.25000
H  -3.84800 17.16600 42.23200
H  -2.93200 17.88100 46.97400
H  -4.58000 17.65400 47.45400
H  -4.05900 19.20000 46.91000


