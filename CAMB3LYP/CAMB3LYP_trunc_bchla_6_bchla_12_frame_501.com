%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 16.95800 -2.04300 28.05800
C  15.67900 0.03700 30.54100
C  18.78700 -3.72100 30.44000
C  18.31300 -3.78900 25.50500
C  15.32700 0.02400 25.64200
N  17.25900 -1.76300 30.19400
C  16.62800 -0.88400 31.01900
C  16.98400 -1.14100 32.48900
C  18.10000 -2.26900 32.35300
C  18.02400 -2.61400 30.86200
C  19.48200 -1.75800 32.74200
C  15.83300 -1.72700 33.38200
C  15.87500 -1.41800 34.85600
H  17.10800 -0.90300 35.45200
N  18.26300 -3.57700 28.00200
C  18.92400 -4.12800 29.11800
C  19.62000 -5.29300 28.63800
C  19.47300 -5.38800 27.23100
C  18.65200 -4.21500 26.84200
C  20.55500 -6.11300 29.50100
C  19.94500 -6.39800 26.27100
O  19.91000 -6.24200 25.06900
C  20.68300 -7.69100 26.77500
N  16.72600 -2.04100 25.82800
C  17.40600 -2.83200 25.05300
C  17.12300 -2.53800 23.54900
C  16.12900 -1.39100 23.56700
C  16.05500 -1.03300 25.09500
C  18.49500 -2.31000 22.75600
C  14.77700 -1.75000 22.91500
C  14.57500 -1.34500 21.40700
N  15.59300 -0.54700 28.03400
C  15.01200 0.20800 26.98500
C  14.13500 1.19500 27.48900
C  14.45500 1.19400 28.91300
C  15.31700 0.08300 29.13900
C  13.41400 2.25000 26.67000
C  14.16900 1.92500 30.09200
O  13.43900 2.88400 30.30900
C  14.95100 1.07800 31.29700
C  15.78000 2.10000 32.01800
O  16.65500 2.79200 31.48700
O  15.36200 2.29400 33.36800
C  16.13600 3.30900 34.14900
H  19.34000 -4.20600 31.24800
H  18.93000 -4.26900 24.74300
H  14.76600 0.65500 24.95000
H  17.39400 -0.20500 32.86900
H  17.90900 -3.13400 32.98800
H  19.41600 -0.81500 33.28400
H  20.16500 -1.70100 31.89400
H  19.97000 -2.35300 33.51400
H  15.82100 -2.80600 33.23200
H  14.89600 -1.40800 32.92500
H  15.72700 -2.36200 35.38000
H  15.08500 -0.70500 35.09200
H  21.59900 -6.03400 29.19800
H  20.20800 -7.14700 29.48400
H  20.60300 -5.83600 30.55400
H  21.71600 -7.39100 26.95300
H  20.76400 -8.41600 25.96500
H  20.25400 -8.04400 27.71200
H  16.67800 -3.47000 23.20000
H  16.55100 -0.48300 23.13900
H  18.49400 -1.34000 22.25900
H  18.63800 -3.00100 21.92500
H  19.26300 -2.34000 23.53000
H  14.01500 -1.24300 23.50800
H  14.64200 -2.82300 23.05200
H  14.22000 -0.31400 21.38800
H  13.95900 -1.97800 20.76800
H  15.61200 -1.30800 21.07500
H  12.51000 1.73500 26.34600
H  13.99600 2.54300 25.79700
H  13.31000 3.23400 27.12900
H  14.20300 0.59900 31.93000
H  16.03100 4.30700 33.72500
H  17.20500 3.14900 34.29200
H  15.57400 3.26700 35.08200
Mg 46.93600 15.62500 28.05800
C  44.93200 15.29600 30.93700
C  49.09000 17.55900 29.83300
C  48.76000 16.00400 25.23500
C  44.45800 14.02100 26.22900
N  46.95800 16.41700 30.08200
C  46.09700 16.05900 31.12300
C  46.62800 16.60100 32.41600
C  47.62900 17.68000 31.93900
C  47.92100 17.18200 30.50000
C  47.12800 19.15800 31.89100
C  47.27700 15.54600 33.38900
C  46.82300 15.65300 34.81000
H  47.91800 15.57200 35.89000
N  48.69300 16.53500 27.65600
C  49.51600 17.23100 28.55700
C  50.69200 17.65500 27.90800
C  50.64200 17.20000 26.50700
C  49.36600 16.54900 26.42000
C  51.74600 18.54700 28.49600
C  51.62200 17.40800 25.27200
O  51.37900 16.91300 24.19700
C  52.94200 18.14800 25.43500
N  46.59000 15.18000 25.99700
C  47.50700 15.46300 25.03600
C  46.95600 15.01700 23.59800
C  45.70400 14.10600 23.96500
C  45.55200 14.44900 25.50800
C  46.57100 16.16300 22.62200
C  45.92300 12.57300 23.70500
C  44.84800 12.01100 22.87300
N  45.05500 14.86000 28.45700
C  44.20200 14.13300 27.62800
C  43.05900 13.68500 28.37600
C  43.29100 14.14000 29.65800
C  44.47900 14.80000 29.68300
C  41.94600 12.80400 27.75100
C  42.74500 14.20900 31.00100
O  41.65000 13.85600 31.38000
C  43.81000 15.00600 31.90600
C  43.16000 16.20600 32.46400
O  42.63000 17.11600 31.77400
O  43.20400 16.12500 33.82500
C  42.27100 17.04800 34.52400
H  49.71300 18.20600 30.45400
H  49.27600 16.09400 24.27700
H  43.65900 13.58300 25.62800
H  45.89800 17.17200 32.99100
H  48.53800 17.72600 32.53900
H  46.41200 19.40800 32.67500
H  46.69600 19.31300 30.90200
H  47.90500 19.92200 31.89600
H  48.35200 15.68900 33.28100
H  47.12000 14.52600 33.03700
H  46.11600 14.82500 34.86400
H  46.37400 16.63200 34.97700
H  52.68800 18.06100 28.75100
H  51.39900 18.92200 29.45800
H  52.04800 19.43100 27.93500
H  53.46800 18.01800 24.49000
H  53.47200 17.83200 26.33400
H  52.75100 19.17700 25.74100
H  47.74100 14.50200 23.04400
H  44.84000 14.47300 23.40900
H  46.63200 17.17800 23.01600
H  45.55400 16.07500 22.23900
H  47.24700 16.08600 21.77000
H  45.87800 12.05000 24.66100
H  46.85000 12.32800 23.18700
H  45.20700 11.93900 21.84600
H  43.86300 12.47600 22.89300
H  44.58300 11.02200 23.24700
H  42.31900 12.37000 26.82300
H  41.00000 13.34300 27.68200
H  41.83700 11.97900 28.45400
H  44.19600 14.38400 32.71400
H  42.26000 16.83600 35.59300
H  41.27900 16.75800 34.17700
H  42.55700 18.06000 34.23600


