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


