%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg 9.20000 47.36600 24.88000
C  7.02600 47.76900 27.74700
C  11.57000 48.44700 26.79100
C  10.94600 47.19300 22.09000
C  6.29100 46.96700 22.73800
N  9.27600 47.92100 27.10800
C  8.32000 48.03400 28.09600
C  8.93500 48.52700 29.35600
C  10.25500 49.18700 28.88100
C  10.38000 48.52600 27.48800
C  10.19000 50.72500 28.75600
C  9.27100 47.32000 30.28400
C  8.91000 47.48600 31.81100
H  8.76300 46.21000 32.73200
N  11.12300 47.75000 24.52700
C  11.99300 48.08000 25.44400
C  13.30700 48.11300 24.83000
C  13.09300 47.62600 23.51900
C  11.67000 47.53200 23.33600
C  14.56100 48.44000 25.66100
C  14.15200 47.33300 22.49300
O  13.84400 46.95300 21.40300
C  15.59100 47.42400 22.88400
N  8.59400 47.38800 22.62900
C  9.59400 47.15300 21.78200
C  9.10700 46.72400 20.44900
C  7.56400 46.85800 20.64800
C  7.41500 47.09900 22.12700
C  9.70100 47.44000 19.18300
C  6.80800 45.57800 20.11700
C  6.00900 45.65200 18.78600
N  7.17400 47.35700 25.08600
C  6.10400 47.11400 24.19600
C  4.83100 47.11900 24.94600
C  5.17200 47.34400 26.27300
C  6.54300 47.52700 26.33600
C  3.42100 47.05000 24.31100
C  4.55900 47.48400 27.63300
O  3.43200 47.38700 28.05800
C  5.73500 47.94200 28.63300
C  5.68600 47.23700 29.86100
O  6.32700 46.19300 30.04600
O  4.84600 47.74900 30.84500
C  4.90400 47.16700 32.22100
H  12.46500 48.75400 27.33500
H  11.53100 47.07200 21.17600
H  5.40800 46.83900 22.11000
H  8.32300 49.32300 29.78100
H  11.13800 49.00400 29.49300
H  10.92100 51.19300 29.41600
H  9.18200 51.02300 29.04300
H  10.33700 51.10600 27.74500
H  10.35100 47.19900 30.19800
H  8.78700 46.45500 29.83100
H  8.00100 48.08400 31.87000
H  9.74500 48.08400 32.17800
H  15.18900 47.54900 25.70000
H  14.38200 48.83100 26.66200
H  15.20000 49.13500 25.11700
H  15.85600 46.71400 23.66700
H  15.91500 48.35300 23.35300
H  16.20400 47.15000 22.02500
H  9.34100 45.66500 20.34000
H  7.22100 47.78300 20.18400
H  10.49600 48.09000 19.54700
H  9.00000 47.94700 18.51900
H  10.12400 46.67900 18.52700
H  6.07000 45.27200 20.85900
H  7.53000 44.76900 20.22800
H  6.23700 44.78300 18.16900
H  6.25500 46.57800 18.26600
H  4.93600 45.62800 18.98000
H  3.49200 46.43400 23.41500
H  3.14200 48.04200 23.95500
H  2.68000 46.64400 25.00000
H  5.61100 49.00700 28.82700
H  5.88500 47.32300 32.67200
H  4.47300 46.16700 32.25100
H  4.32200 47.83100 32.86000


