%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 44.01800 2.91800 47.66100
C  42.06000 5.71500 47.01900
C  41.00600 1.01600 47.31200
C  45.79900 0.15700 47.76200
C  46.87900 4.91700 47.24100
N  41.83200 3.33000 47.09800
C  41.29400 4.59200 46.98100
C  39.79500 4.50000 46.68300
C  39.48700 3.02600 47.06600
C  40.84700 2.38200 47.12900
C  38.69200 2.74400 48.40100
C  39.44800 4.83900 45.22400
C  38.00300 4.93400 44.91400
H  37.59500 5.76800 43.64000
N  43.49300 0.83400 47.54100
C  42.21400 0.26700 47.42800
C  42.35900 -1.16200 47.43400
C  43.77600 -1.43800 47.42700
C  44.43200 -0.09700 47.60800
C  41.15600 -2.13400 47.28600
C  44.34200 -2.85900 47.34100
O  43.62700 -3.84400 47.42200
C  45.80100 -3.04100 47.16300
N  46.06200 2.62900 47.74100
C  46.53800 1.38200 47.82400
C  48.11500 1.40300 47.69600
C  48.51000 2.91300 47.47000
C  47.02000 3.56100 47.42600
C  48.93400 0.68700 48.80200
C  49.24600 3.07700 46.12300
C  50.25900 4.22700 46.06800
N  44.47400 4.92000 47.22100
C  45.65200 5.61300 47.20700
C  45.34700 7.04400 47.09000
C  43.98000 7.11700 47.10900
C  43.47900 5.83500 47.12000
C  46.36600 8.18800 47.10600
C  42.87400 8.02300 47.02800
O  42.81600 9.22700 46.97900
C  41.57700 7.14800 46.98100
C  40.79100 7.44500 48.21000
O  41.19100 7.46200 49.34100
O  39.43700 7.57100 47.85400
C  38.52600 7.78100 48.94300
H  40.02600 0.53600 47.28500
H  46.51600 -0.66500 47.82200
H  47.74600 5.57800 47.17700
H  39.30800 5.14400 47.41400
H  38.96600 2.43700 46.31200
H  37.88000 2.08800 48.08800
H  38.53600 3.57700 49.08700
H  39.33000 2.09900 49.00500
H  39.75400 4.10000 44.48300
H  39.93900 5.76800 44.93400
H  37.52300 5.35800 45.79600
H  37.59600 3.92700 44.82400
H  41.34100 -2.41600 46.24900
H  40.16900 -1.67200 47.25600
H  41.01800 -3.02000 47.90600
H  46.05900 -2.38200 46.33400
H  45.98200 -4.06600 46.83900
H  46.21200 -2.73700 48.12600
H  48.33800 0.83400 46.79300
H  49.07900 3.37000 48.27900
H  49.32800 -0.23900 48.38200
H  48.30000 0.56800 49.68000
H  49.82500 1.27000 49.03500
H  48.45300 3.22300 45.39000
H  49.73200 2.16800 45.76800
H  51.18600 3.72000 45.79800
H  50.32400 4.81900 46.98100
H  50.02100 4.90100 45.24500
H  46.22000 8.97700 46.36900
H  47.30400 7.65100 46.96300
H  46.44100 8.50600 48.14600
H  41.03900 7.44900 46.08200
H  38.18200 8.81200 48.86000
H  38.86300 7.60600 49.96500
H  37.69400 7.09100 48.80500
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


