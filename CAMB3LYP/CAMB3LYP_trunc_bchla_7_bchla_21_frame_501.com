%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 25.26700 0.12000 29.48500
C  27.14100 -0.08700 32.36300
C  22.44200 0.65100 31.48400
C  23.33900 0.24400 26.71700
C  27.97000 -0.76800 27.63500
N  24.82700 -0.00300 31.72000
C  25.83000 0.13900 32.68700
C  25.26100 0.57400 34.05200
C  23.73000 0.68200 33.70700
C  23.61100 0.48500 32.21100
C  22.87000 -0.30200 34.49000
C  25.80700 1.87500 34.74000
C  26.64900 1.68700 36.11600
H  25.98700 2.43400 37.27800
N  23.21300 0.29800 29.14100
C  22.23900 0.64000 30.09300
C  20.97400 1.01400 29.45200
C  21.16000 0.82900 27.97700
C  22.63600 0.50000 27.85200
C  19.74900 1.44500 30.22800
C  20.23300 0.96600 26.86300
O  20.59800 0.66100 25.71800
C  18.79000 1.40700 27.03500
N  25.61300 -0.29200 27.55300
C  24.69700 -0.04800 26.51300
C  25.29200 -0.18000 25.05300
C  26.81800 -0.54700 25.43400
C  26.79800 -0.51500 26.97900
C  24.57500 -1.29000 24.13700
C  27.83300 0.46900 24.80600
C  29.13600 -0.06400 24.29800
N  27.18000 -0.44400 29.94400
C  28.21700 -0.73500 29.05300
C  29.45900 -0.85000 29.73800
C  29.08500 -0.49600 31.08700
C  27.69300 -0.34600 31.12700
C  30.82000 -1.11200 29.24500
C  29.58400 -0.41400 32.42400
O  30.70000 -0.50100 32.94800
C  28.28000 -0.28500 33.33500
C  28.18400 -1.48700 34.20500
O  27.52600 -2.44100 34.06000
O  28.92900 -1.30000 35.36700
C  28.48600 -2.07100 36.51000
H  21.58300 0.90400 32.10900
H  22.84400 0.32400 25.74700
H  28.87000 -0.90500 27.03200
H  25.38700 -0.24700 34.75700
H  23.31500 1.64400 34.00900
H  22.43700 -1.01900 33.79300
H  22.09100 0.33900 34.90100
H  23.35600 -0.79100 35.33400
H  24.95600 2.52500 34.94500
H  26.49900 2.40700 34.08700
H  27.58200 2.20000 35.88200
H  26.76100 0.62900 36.35200
H  19.34600 0.54600 30.69600
H  19.02800 1.84900 29.51900
H  19.97500 2.10600 31.06500
H  18.28300 1.10300 26.11900
H  18.67500 2.47500 27.22300
H  18.18400 0.86400 27.76000
H  25.32400 0.80000 24.57800
H  27.16700 -1.50900 25.05800
H  25.34100 -1.79600 23.54900
H  23.91500 -0.92500 23.35100
H  24.00900 -2.02000 24.71600
H  28.05800 1.32200 25.44600
H  27.45300 0.92700 23.89300
H  29.88300 0.35500 24.97200
H  29.54600 0.28700 23.35100
H  29.13600 -1.15300 24.35500
H  30.88400 -0.76300 28.21400
H  30.89900 -2.19900 29.25700
H  31.53400 -0.73600 29.97800
H  28.35300 0.59900 33.96800
H  28.85600 -3.07500 36.30500
H  27.39700 -2.05300 36.47000
H  28.92200 -1.70000 37.43800
Mg 15.62100 51.82400 26.04800
C  17.47400 50.60700 28.62800
C  13.11000 52.75700 28.20700
C  13.81500 52.82700 23.36600
C  17.99400 50.42600 23.83300
N  15.23600 51.49000 28.16800
C  16.26900 51.27600 29.02600
C  15.94600 51.61800 30.44100
C  14.40600 52.07300 30.26100
C  14.21100 52.11900 28.78100
C  13.38100 51.07100 30.96200
C  16.88200 52.82900 30.90600
C  17.11300 53.00000 32.40900
H  16.27700 52.29800 33.50800
N  13.77400 52.76600 25.77300
C  12.87200 53.04400 26.81300
C  11.70500 53.64500 26.21000
C  11.90800 53.69900 24.81600
C  13.23300 53.12400 24.62900
C  10.53200 54.15900 27.07700
C  10.86200 54.19000 23.68400
O  11.18100 54.33400 22.49700
C  9.41300 54.55600 23.96000
N  15.93700 51.69300 23.86100
C  14.99700 52.16200 23.00600
C  15.27300 51.71300 21.57200
C  16.82600 51.21000 21.68000
C  17.02200 51.18500 23.20400
C  14.21700 50.69700 20.93700
C  17.94900 52.21900 21.11400
C  19.01900 51.70400 20.16900
N  17.42400 50.77200 26.12100
C  18.25400 50.26800 25.19100
C  19.45700 49.67400 25.85400
C  19.17000 49.86100 27.22700
C  17.86900 50.41600 27.30000
C  20.69200 49.07400 25.28100
C  19.74000 49.61300 28.49800
O  20.83000 49.08800 28.80500
C  18.66300 50.08200 29.47300
C  18.25000 48.96100 30.28600
O  17.83600 47.88700 29.97900
O  18.30900 49.39800 31.60600
C  17.88900 48.52200 32.76600
H  12.38500 53.27700 28.83600
H  13.18500 53.06500 22.50600
H  18.65800 49.80900 23.22500
H  16.01600 50.87800 31.23800
H  14.23400 53.06400 30.68100
H  12.50100 51.52900 31.41300
H  13.92900 50.42000 31.64300
H  12.99700 50.42600 30.17100
H  16.51300 53.78500 30.53400
H  17.84200 52.64500 30.42300
H  17.07200 54.07500 32.58900
H  18.14900 52.74600 32.63000
H  10.69700 54.07400 28.15100
H  9.65000 53.59900 26.76800
H  10.36400 55.20000 26.80000
H  9.47200 55.48000 24.53400
H  8.95500 53.91200 24.71100
H  8.82100 54.48200 23.04800
H  15.27500 52.60400 20.94300
H  16.96000 50.24900 21.18300
H  14.69800 49.77000 20.62400
H  13.70900 51.16400 20.09400
H  13.48400 50.38900 21.68300
H  18.44000 52.77000 21.91600
H  17.40700 53.00700 20.59200
H  19.97700 51.93800 20.63200
H  18.94100 52.12500 19.16700
H  18.97400 50.61500 20.12300
H  20.76700 49.12700 24.19500
H  20.76600 48.01700 25.53800
H  21.55200 49.59200 25.70500
H  19.01800 50.87900 30.12700
H  17.50200 47.57000 32.40200
H  17.07000 48.99900 33.30600
H  18.78800 48.39000 33.36800


