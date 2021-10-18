%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 47.67100 15.89000 28.45400
C  45.34400 15.65500 31.04800
C  49.48000 18.07500 30.51500
C  49.68700 16.65900 25.89600
C  45.30600 14.62100 26.30000
N  47.50700 16.71700 30.64000
C  46.40900 16.47700 31.40500
C  46.52000 17.23500 32.78700
C  47.75800 18.22000 32.48000
C  48.27100 17.69900 31.12700
C  47.16600 19.65600 32.25100
C  46.90400 16.16800 33.87700
C  47.20400 16.68600 35.31000
H  48.31600 16.29900 36.23900
N  49.43000 17.01400 28.32000
C  50.06900 17.82700 29.27100
C  51.20000 18.39100 28.66300
C  51.38300 17.85100 27.39700
C  50.14500 17.12800 27.15700
C  52.11500 19.27100 29.45000
C  52.51500 18.07600 26.34100
O  52.57400 17.52300 25.24800
C  53.78800 18.86100 26.66900
N  47.44200 15.76700 26.39600
C  48.48400 16.01700 25.57200
C  48.26000 15.53000 24.12700
C  46.80100 15.00400 24.22400
C  46.51000 15.00800 25.71400
C  48.61100 16.58300 23.02400
C  46.48600 13.65000 23.53100
C  46.96600 12.43800 24.21000
N  45.75600 15.16700 28.59900
C  44.93800 14.69700 27.68600
C  43.69800 14.18000 28.27400
C  43.86300 14.56100 29.64400
C  45.08900 15.18100 29.78600
C  42.77800 13.18000 27.74800
C  43.21100 14.49100 30.94000
O  42.19200 13.93500 31.28200
C  44.14500 15.27600 31.86900
C  43.41900 16.43800 32.55500
O  43.17000 17.49000 32.02600
O  43.04200 16.05500 33.82000
C  42.33900 17.13300 34.55300
H  50.09800 18.64700 31.20900
H  50.41400 16.70200 25.08200
H  44.81200 14.14000 25.45400
H  45.63400 17.83000 33.01000
H  48.52200 18.24600 33.25800
H  47.13800 20.01300 31.22200
H  47.84000 20.42900 32.62100
H  46.16600 19.79300 32.66200
H  47.75800 15.61000 33.49400
H  46.09800 15.43400 33.89800
H  46.24800 16.58700 35.82500
H  47.42100 17.73200 35.09700
H  51.80300 19.48300 30.47300
H  52.34700 20.14300 28.83700
H  53.02500 18.67700 29.53900
H  53.61800 19.92500 26.83400
H  54.61200 18.70300 25.97300
H  54.19200 18.41700 27.57900
H  49.00300 14.73300 24.13500
H  46.20400 15.80500 23.78800
H  49.40000 17.25400 23.36200
H  47.72600 17.18100 22.80500
H  48.98600 16.12600 22.10800
H  46.99700 13.77600 22.57700
H  45.42300 13.59700 23.29500
H  46.25100 12.14600 24.97900
H  47.97200 12.58100 24.60400
H  46.98800 11.57800 23.54000
H  43.37600 12.59600 27.04800
H  41.93300 13.49000 27.13300
H  42.41600 12.61300 28.60500
H  44.52600 14.58100 32.61800
H  41.57800 16.67700 35.18600
H  41.70900 17.79000 33.95300
H  43.02300 17.66800 35.21100
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


