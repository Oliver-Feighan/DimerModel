%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
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


