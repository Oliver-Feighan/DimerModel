%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg -1.96800 16.99100 27.07800
C  -2.18500 14.77300 29.87900
C  -3.22300 19.43300 29.10400
C  -2.21300 18.84100 24.51400
C  -1.82200 14.16700 25.04100
N  -2.58500 16.98200 29.20400
C  -2.56200 16.08300 30.18300
C  -2.90400 16.68200 31.58000
C  -3.74100 17.95500 31.11000
C  -3.12300 18.15900 29.73800
C  -5.23100 17.67600 30.94100
C  -1.62300 17.10100 32.37700
C  -1.84400 17.38800 33.86500
H  -0.80000 17.02400 34.81100
N  -2.15300 18.95200 26.92600
C  -2.67100 19.79100 27.85500
C  -2.58800 21.14700 27.27100
C  -2.15500 21.00600 25.85300
C  -2.18200 19.57400 25.69800
C  -2.85500 22.39700 28.05200
C  -1.78700 22.00400 24.80100
O  -1.66500 21.68100 23.57000
C  -1.65800 23.40500 25.24100
N  -2.23300 16.57800 25.04900
C  -2.24700 17.56000 24.14100
C  -2.16400 17.01500 22.74600
C  -1.84400 15.45800 22.89800
C  -1.88200 15.35900 24.44300
C  -3.41900 17.19600 21.88100
C  -0.59500 14.92100 22.15500
C  0.65800 15.50100 22.77200
N  -1.81800 14.84700 27.38100
C  -1.79500 13.85500 26.39000
C  -1.79500 12.55300 27.08200
C  -1.83500 12.88000 28.43600
C  -1.88500 14.29600 28.58000
C  -1.81900 11.27600 26.45500
C  -1.87200 12.36200 29.79400
O  -1.76000 11.24400 30.29900
C  -2.22300 13.56300 30.77900
C  -1.25100 13.50900 31.92600
O  -0.03700 13.73600 31.74700
O  -1.93000 13.43500 33.08100
C  -1.02500 13.40800 34.26100
H  -3.82700 20.21100 29.57500
H  -2.18600 19.44400 23.60400
H  -1.69200 13.24400 24.47300
H  -3.52200 16.04200 32.21000
H  -3.55300 18.82200 31.74300
H  -5.45200 16.61400 31.04700
H  -5.47900 18.01800 29.93600
H  -5.82000 18.29200 31.62100
H  -1.29900 18.08700 32.04400
H  -0.83800 16.35200 32.26700
H  -2.81200 17.01700 34.20300
H  -1.90600 18.46800 33.99900
H  -3.31200 22.17700 29.01600
H  -3.61900 22.95400 27.51000
H  -1.96100 23.01100 28.16200
H  -2.60300 23.73500 25.67300
H  -1.30300 23.96100 24.37400
H  -0.92000 23.45700 26.04200
H  -1.36000 17.49200 22.18700
H  -2.66000 14.83500 22.53100
H  -3.23500 17.73100 20.94900
H  -4.17200 17.70500 22.48200
H  -3.93500 16.27100 21.62500
H  -0.61200 15.23200 21.11100
H  -0.51600 13.83400 22.17600
H  0.74600 16.46900 22.28000
H  1.51900 14.87000 22.54900
H  0.52700 15.68800 23.83800
H  -1.10100 10.65800 26.99400
H  -1.56400 11.34900 25.39800
H  -2.85900 10.99200 26.61500
H  -3.25100 13.41000 31.10600
H  -1.27000 12.51300 34.83300
H  -1.13500 14.37200 34.75800
H  0.03200 13.23700 34.05400
Mg 29.46400 59.14500 41.03800
C  26.57800 57.38900 40.28900
C  31.34900 56.53100 39.78200
C  32.15000 61.18900 41.08600
C  27.35400 61.82000 41.94900
N  28.99200 57.29200 39.87100
C  27.77200 56.74800 39.89300
C  27.77900 55.38900 39.12700
C  29.32900 54.99700 39.21800
C  29.96700 56.34400 39.73200
C  29.58100 53.72700 40.06600
C  27.37400 55.59400 37.68500
C  26.46000 54.52700 37.10300
H  26.57800 54.02100 35.61600
N  31.53300 58.78800 40.65000
C  32.11000 57.71400 40.08400
C  33.51500 57.95800 40.03300
C  33.80600 59.24700 40.57400
C  32.47900 59.81600 40.84000
C  34.47700 56.89100 39.66100
C  35.16200 59.85800 40.59300
O  36.13000 59.21400 40.22400
C  35.39700 61.26400 41.13200
N  29.71300 61.28200 41.56500
C  30.94400 61.84600 41.41400
C  30.85900 63.33400 41.52000
C  29.29200 63.52200 41.79300
C  28.69900 62.13600 41.84900
C  31.79900 63.96200 42.56600
C  28.51800 64.48000 40.78300
C  27.64200 65.63900 41.31400
N  27.43700 59.54700 41.14000
C  26.74400 60.62700 41.58900
C  25.41200 60.28000 41.81900
C  25.25400 58.98400 41.34800
C  26.48900 58.61100 40.88000
C  24.32700 61.25500 42.36300
C  24.37200 57.92800 41.16800
O  23.15900 57.91800 41.34600
C  25.14900 56.75500 40.39700
C  25.25800 55.58200 41.18700
O  25.99300 55.31600 42.14100
O  24.33400 54.68200 40.66500
C  24.17500 53.35900 41.24400
H  31.95200 55.65500 39.53500
H  32.91900 61.96200 41.15300
H  26.80500 62.68200 42.33400
H  27.08300 54.71900 39.63100
H  29.81700 54.88800 38.25000
H  28.66900 53.37400 40.54700
H  30.34000 53.93800 40.82000
H  30.00600 52.97200 39.40500
H  28.16800 55.74800 36.95400
H  26.92700 56.58800 37.67000
H  25.44000 54.85100 37.30900
H  26.66400 53.69100 37.77200
H  34.81300 56.52500 40.63100
H  35.29700 57.28300 39.05900
H  34.02200 56.08300 39.08800
H  35.10600 61.93900 40.32700
H  36.47900 61.37300 41.20200
H  34.84400 61.42100 42.05800
H  31.14300 63.61000 40.50500
H  29.10900 64.04800 42.73000
H  32.59600 64.44300 41.99800
H  32.19100 63.23400 43.27700
H  31.26900 64.75100 43.09900
H  27.85100 63.86700 40.17700
H  29.18800 65.01800 40.11200
H  26.69000 65.11600 41.40400
H  27.54400 66.41700 40.55800
H  27.77900 66.00900 42.33000
H  24.06200 61.95200 41.56900
H  24.63400 61.90700 43.18100
H  23.42900 60.71200 42.65900
H  24.61600 56.61800 39.45600
H  23.13300 53.07400 41.10200
H  24.52900 53.25500 42.27000
H  24.72400 52.65900 40.61400


