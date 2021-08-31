%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 44.80500 2.87200 46.67800
C  42.56400 5.57100 46.81800
C  42.13800 0.79200 46.35700
C  47.01300 0.52400 46.12300
C  47.34200 5.37800 46.17800
N  42.64100 3.17100 46.66900
C  41.94200 4.34400 46.71700
C  40.48300 4.04900 46.33800
C  40.34000 2.59000 46.75700
C  41.78900 2.14300 46.48400
C  39.96100 2.42600 48.24700
C  40.23600 4.25700 44.79800
C  39.26900 5.45400 44.33400
H  39.83900 6.57400 43.44700
N  44.60700 0.93700 46.23600
C  43.41900 0.24600 46.16900
C  43.72000 -1.14400 46.02500
C  45.14200 -1.32500 46.06900
C  45.67600 0.06600 46.14000
C  42.62100 -2.25300 45.89000
C  45.87200 -2.71700 46.06400
O  45.23700 -3.72500 46.11700
C  47.31400 -2.74100 46.27700
N  46.83800 2.95400 45.90000
C  47.54500 1.82700 46.15300
C  49.01900 2.09000 46.41200
C  49.17500 3.56200 46.06200
C  47.67600 4.02200 46.03300
C  49.44600 1.78300 47.89500
C  49.96600 3.97900 44.73400
C  51.33300 4.60700 44.78400
N  45.00300 4.99200 46.72600
C  46.07700 5.86500 46.51000
C  45.65700 7.24900 46.55200
C  44.29400 7.16900 46.69900
C  43.89900 5.79300 46.77900
C  46.53700 8.41600 46.37700
C  43.08200 7.96600 46.81100
O  42.95000 9.19400 46.88400
C  41.94500 6.98400 47.00000
C  41.32800 7.26400 48.36600
O  42.00600 7.20100 49.38900
O  40.00100 7.70700 48.25300
C  39.36700 8.02300 49.50500
H  41.26200 0.14300 46.41400
H  47.71200 -0.30200 46.26300
H  48.13300 6.13000 46.13700
H  39.87300 4.80700 46.83000
H  39.59400 2.12600 46.11300
H  39.12500 1.73800 48.12600
H  39.66800 3.42600 48.56800
H  40.73100 2.03300 48.91100
H  39.76200 3.35000 44.42300
H  41.21000 4.30800 44.31200
H  38.98200 5.96200 45.25500
H  38.33700 5.13200 43.87100
H  42.71600 -2.97100 45.07600
H  41.68000 -1.71700 45.76200
H  42.51800 -2.85500 46.79300
H  47.88700 -2.17400 45.54300
H  47.64600 -3.77900 46.25800
H  47.54400 -2.49600 47.31300
H  49.67000 1.44400 45.82200
H  49.62700 4.10200 46.89400
H  49.91800 0.80800 48.01500
H  48.60400 1.96800 48.56200
H  50.20200 2.52100 48.16700
H  49.41200 4.48800 43.94600
H  50.09500 2.95900 44.37100
H  52.11800 3.97300 44.37300
H  51.60000 4.81700 45.82000
H  51.22600 5.52900 44.21300
H  47.39900 8.09600 45.79100
H  46.88900 8.66700 47.37800
H  45.99500 9.26600 45.96500
H  41.25100 7.27600 46.21300
H  39.04600 7.15600 50.08300
H  38.48400 8.64500 49.36100
H  39.98200 8.61300 50.18400
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


