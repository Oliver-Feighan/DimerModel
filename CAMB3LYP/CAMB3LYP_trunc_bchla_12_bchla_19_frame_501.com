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
Mg 25.47400 50.80800 26.36000
C  23.42800 51.77500 29.14000
C  27.78000 49.80000 28.46900
C  27.48600 50.74800 23.72700
C  22.89100 52.35600 24.24800
N  25.55500 50.75400 28.55700
C  24.64000 51.22200 29.52100
C  25.12200 51.01800 31.00700
C  26.56400 50.28500 30.68300
C  26.67300 50.24400 29.16200
C  27.79400 51.06900 31.31300
C  24.14500 50.22000 31.93900
C  24.75600 49.66700 33.22100
H  24.21600 50.30100 34.51600
N  27.35300 50.24000 26.12300
C  28.22400 49.80700 27.09500
C  29.50900 49.42000 26.55200
C  29.38500 49.70200 25.13400
C  27.97500 50.19600 24.94200
C  30.61900 48.81900 27.27800
C  30.38700 49.50900 24.00900
O  30.13300 49.66800 22.84800
C  31.78500 49.18000 24.46100
N  25.24600 51.52100 24.28400
C  26.25500 51.32600 23.39800
C  25.83300 51.76400 21.97000
C  24.32300 52.13900 22.16700
C  24.10500 51.95100 23.65600
C  26.64800 52.90800 21.27600
C  23.28700 51.29000 21.34900
C  23.20200 51.69600 19.86100
N  23.59000 51.86200 26.56100
C  22.66400 52.33500 25.66000
C  21.52300 52.77100 26.36200
C  21.75600 52.64400 27.68000
C  23.06700 52.09200 27.80400
C  20.26700 53.21900 25.62600
C  21.30400 52.86100 29.05800
O  20.28500 53.34300 29.48100
C  22.42600 52.30600 30.05900
C  23.02400 53.48800 30.75300
O  23.93000 54.23000 30.44700
O  22.33300 53.57900 31.94500
C  22.72800 54.66700 32.87100
H  28.62400 49.48900 29.08800
H  28.23800 50.84300 22.94000
H  22.18300 52.67600 23.48000
H  25.36700 52.00000 31.41200
H  26.54900 49.29900 31.14800
H  28.47500 50.35700 31.77900
H  27.60200 51.88000 32.01500
H  28.41300 51.49500 30.52400
H  23.78800 49.47000 31.23200
H  23.21100 50.69500 32.24000
H  25.83400 49.82400 33.19700
H  24.45000 48.62700 33.33100
H  30.63000 47.83900 26.80100
H  30.42700 48.80700 28.35100
H  31.61000 49.25400 27.15000
H  31.91600 48.09800 24.46700
H  32.15000 49.82600 25.26000
H  32.40300 49.70200 23.73100
H  25.75100 50.83500 21.40500
H  24.20500 53.21400 22.02900
H  26.08400 53.79400 20.98300
H  26.95200 52.49600 20.31400
H  27.50400 53.28000 21.83900
H  22.29300 51.43100 21.77200
H  23.59200 50.25400 21.50300
H  23.57200 52.71200 19.72200
H  22.15600 51.57400 19.58000
H  23.82700 51.09700 19.19800
H  19.37700 53.06800 26.23700
H  20.20500 52.68100 24.68000
H  20.36300 54.25800 25.31000
H  22.01700 51.57800 30.76000
H  21.79500 55.04700 33.28700
H  23.27800 55.48100 32.39900
H  23.31400 54.08500 33.58300


