%nproc=24
%mem=175gb
#p wB97XD/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 1.29600 7.71300 26.88100
C  1.46900 9.82700 29.67000
C  2.42600 5.06900 28.94700
C  1.47300 5.72500 24.13400
C  0.92500 10.50000 24.87700
N  1.69900 7.49100 29.04900
C  1.51500 8.49600 30.00300
C  1.43800 7.84300 31.43000
C  1.83300 6.38800 31.15200
C  2.05200 6.30500 29.60100
C  3.12000 6.01700 31.93700
C  0.04000 8.06800 31.99300
C  -0.12000 8.13500 33.50300
H  1.05800 7.66300 34.38500
N  1.84800 5.68100 26.53500
C  2.31100 4.77600 27.50900
C  2.63500 3.49600 26.90400
C  2.31300 3.59100 25.49300
C  1.87500 5.04400 25.33500
C  3.03900 2.24400 27.69300
C  2.37500 2.55400 24.38800
O  2.03400 2.84400 23.24900
C  2.64700 1.14500 24.75200
N  1.13700 8.09300 24.78900
C  1.15400 7.08500 23.87300
C  0.77800 7.64200 22.47500
C  0.53200 9.18300 22.70900
C  0.76500 9.29400 24.23800
C  1.89300 7.33600 21.45900
C  -0.89600 9.67800 22.35700
C  -0.86300 10.38200 21.00400
N  1.30900 9.79300 27.10700
C  1.22300 10.82400 26.23200
C  1.50800 12.06900 26.98600
C  1.50300 11.75700 28.31400
C  1.42300 10.34100 28.32600
C  1.55300 13.42200 26.29600
C  1.55800 12.26300 29.69200
O  1.68500 13.41500 30.08900
C  1.44800 11.03500 30.57500
C  2.56000 11.13900 31.53300
O  3.76600 10.89700 31.38900
O  2.04800 11.65700 32.70500
C  2.99200 11.90900 33.77200
H  2.69600 4.26200 29.63100
H  1.45200 4.99400 23.32300
H  0.89300 11.33000 24.16800
H  2.21800 8.40600 31.94400
H  1.06000 5.77700 31.61800
H  4.00400 5.85500 31.32000
H  2.97100 4.99800 32.29300
H  3.38200 6.77100 32.67900
H  -0.56100 7.28300 31.53300
H  -0.37400 8.97400 31.55200
H  -1.02500 7.59600 33.78000
H  -0.31800 9.19800 33.64500
H  3.22100 2.42800 28.75200
H  3.95300 1.82100 27.27400
H  2.28900 1.47200 27.52500
H  1.81900 0.77800 25.35900
H  3.66100 1.02700 25.13400
H  2.61600 0.48900 23.88200
H  -0.10600 7.10300 22.13500
H  1.34500 9.68200 22.18000
H  2.33200 8.19000 20.94400
H  1.31000 6.77100 20.73200
H  2.77000 6.84600 21.88100
H  -1.17200 10.51700 22.99600
H  -1.64500 8.88600 22.34300
H  -1.71200 9.95300 20.47200
H  -0.02300 10.09700 20.37100
H  -0.94100 11.46900 21.03900
H  1.59500 13.38900 25.20700
H  2.38500 14.03800 26.63500
H  0.58900 13.91000 26.44400
H  0.47700 11.06900 31.06800
H  3.81000 12.49700 33.35500
H  3.40600 11.00100 34.21100
H  2.51900 12.68100 34.37800
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


