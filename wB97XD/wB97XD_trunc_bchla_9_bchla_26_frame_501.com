%nproc=24
%mem=175gb
#p wB97XD/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 35.26200 1.75200 29.96200
C  32.85500 2.58300 32.32000
C  37.64200 1.80100 32.35100
C  37.51600 1.12500 27.56800
C  32.74800 1.88400 27.52800
N  35.16600 2.13300 32.11400
C  34.11900 2.35200 32.87300
C  34.63200 2.45300 34.38300
C  36.03300 2.95200 34.09900
C  36.35200 2.21400 32.76000
C  36.05200 4.49500 34.04900
C  34.37500 1.08800 35.15500
C  34.37600 1.20900 36.68200
H  34.84400 2.50900 37.31300
N  37.24000 1.32700 29.98300
C  38.04600 1.32500 31.11900
C  39.36800 0.74800 30.73100
C  39.37200 0.59900 29.27900
C  38.02100 1.10700 28.92100
C  40.46600 0.47300 31.73100
C  40.47100 0.00700 28.28100
O  40.29400 -0.12300 27.06500
C  41.78300 -0.42500 28.90200
N  35.17700 1.65000 27.89700
C  36.22800 1.39600 27.12200
C  35.89000 1.47100 25.60100
C  34.35600 1.46600 25.66700
C  34.03700 1.61900 27.12100
C  36.49900 2.69400 24.86000
C  33.56700 0.32100 24.88900
C  33.00100 -0.79800 25.70000
N  33.26100 2.34500 29.84600
C  32.36300 2.28700 28.83400
C  31.06300 2.71100 29.21500
C  31.19600 2.84400 30.63300
C  32.52400 2.60100 30.90000
C  29.78600 2.98300 28.46300
C  30.54600 3.13100 31.87700
O  29.40100 3.42500 32.14500
C  31.52900 2.83000 33.05300
C  31.55800 3.97100 33.97300
O  31.92600 5.08300 33.62600
O  30.95500 3.71000 35.23000
C  30.85200 4.85600 36.14700
H  38.39100 1.89200 33.13900
H  38.28200 0.96200 26.80800
H  32.00200 1.96400 26.73500
H  34.02600 3.17400 34.93000
H  36.78800 2.59500 34.79900
H  36.64600 4.74400 34.92900
H  35.05700 4.93600 34.12200
H  36.63900 4.75700 33.16900
H  35.14200 0.35900 34.89200
H  33.36400 0.79400 34.87200
H  34.98900 0.39000 37.05700
H  33.35800 1.00500 37.01300
H  40.05900 0.17200 32.69700
H  41.18400 1.27200 31.91300
H  41.03100 -0.39400 31.38900
H  42.08900 0.44800 29.47800
H  42.58200 -0.61900 28.18700
H  41.63400 -1.33400 29.48500
H  36.25600 0.59300 25.06900
H  34.05000 2.44500 25.30000
H  35.69900 3.06000 24.21600
H  37.40400 2.52000 24.28000
H  36.80800 3.41800 25.61400
H  34.19600 -0.03500 24.07300
H  32.66800 0.75100 24.44800
H  31.91500 -0.71000 25.66100
H  33.14600 -0.73200 26.77900
H  33.24300 -1.79400 25.32900
H  29.91400 2.90200 27.38300
H  29.37400 3.93500 28.79900
H  28.98400 2.25500 28.58700
H  31.23100 1.94400 33.61300
H  30.55000 4.60700 37.16400
H  30.10200 5.48600 35.66800
H  31.81400 5.36300 36.08200
Mg -9.66500 18.33400 42.78400
C  -6.24000 17.80500 42.49600
C  -9.09900 21.67600 42.05300
C  -12.91200 18.87400 42.48300
C  -10.14800 14.91400 43.00500
N  -7.91600 19.56100 42.17100
C  -6.59800 19.16200 42.25800
C  -5.63300 20.29600 42.07600
C  -6.57300 21.55800 42.31900
C  -7.95100 20.93600 42.13600
C  -6.37100 22.54100 43.56800
C  -4.81400 20.28500 40.70300
C  -5.51800 19.69700 39.44400
H  -5.27900 20.48100 38.17600
N  -10.85700 20.01500 42.21400
C  -10.42000 21.32200 41.96500
C  -11.51000 22.18000 41.58900
C  -12.66600 21.36900 41.68400
C  -12.20600 20.01400 42.11100
C  -11.29500 23.68900 41.46200
C  -14.07900 21.93000 41.40500
O  -14.25800 23.08100 40.99400
C  -15.20100 21.03000 41.71200
N  -11.24800 17.01300 42.63400
C  -12.52300 17.52900 42.79400
C  -13.59000 16.43700 43.02800
C  -12.64200 15.17100 42.94600
C  -11.28500 15.70200 42.80600
C  -14.41600 16.47100 44.27500
C  -12.97300 14.21900 41.68800
C  -13.02000 12.77100 42.10100
N  -8.45700 16.64500 42.95800
C  -8.80300 15.30900 43.04100
C  -7.58100 14.47400 43.08900
C  -6.53000 15.43300 42.87800
C  -7.12900 16.69600 42.80500
C  -7.42000 13.00900 43.03700
C  -5.08700 15.65900 42.64900
O  -4.19400 14.85200 42.53200
C  -4.84400 17.18800 42.43700
C  -3.86600 17.68600 43.40500
O  -2.75700 18.18900 43.22300
O  -4.28400 17.22800 44.68200
C  -3.45300 17.67600 45.82000
H  -8.86100 22.74200 42.07500
H  -13.99400 18.96100 42.60100
H  -10.40600 13.86400 43.15400
H  -4.85200 20.25200 42.83500
H  -6.35600 22.13700 41.42100
H  -5.49500 22.34800 44.18700
H  -7.18200 22.50000 44.29500
H  -6.14100 23.54200 43.20200
H  -3.84000 19.81700 40.84200
H  -4.55300 21.32400 40.50200
H  -6.59000 19.59400 39.61500
H  -5.14100 18.69400 39.24600
H  -11.94900 24.16700 40.73200
H  -10.29700 23.87500 41.06600
H  -11.49900 24.10000 42.45100
H  -15.27200 20.92300 42.79400
H  -15.05600 20.07500 41.20600
H  -16.08400 21.54300 41.33100
H  -14.18800 16.36800 42.11900
H  -12.69300 14.77700 43.96100
H  -14.92900 15.51000 44.31600
H  -15.20300 17.20800 44.11500
H  -13.79200 16.69500 45.14100
H  -12.21000 14.42700 40.93800
H  -13.86700 14.45400 41.10900
H  -12.15800 12.27100 41.66000
H  -13.81800 12.22200 41.60000
H  -12.89800 12.69100 43.18100
H  -8.34400 12.43000 43.00100
H  -6.93100 12.73100 43.97000
H  -6.96900 12.89900 42.05100
H  -4.46800 17.26400 41.41700
H  -4.21000 17.98200 46.54300
H  -2.74000 18.49300 45.70900
H  -3.06600 16.75200 46.25000


