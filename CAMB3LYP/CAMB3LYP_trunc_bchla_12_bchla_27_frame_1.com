%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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
Mg -5.21000 24.54100 27.19500
C  -3.54100 26.38400 29.69300
C  -6.38000 22.49900 29.72600
C  -6.78300 22.59100 24.87600
C  -3.73200 26.39300 24.81700
N  -5.00900 24.47500 29.48500
C  -4.39300 25.40300 30.31600
C  -4.77800 25.15700 31.76100
C  -5.38000 23.70500 31.74000
C  -5.61200 23.56800 30.25200
C  -4.48700 22.58300 32.40100
C  -5.82800 26.14800 32.28500
C  -5.66000 26.57200 33.74800
H  -5.92800 25.52000 34.86000
N  -6.48000 22.85300 27.28100
C  -6.76500 22.19100 28.42600
C  -7.52200 21.02600 28.10100
C  -7.81700 21.04500 26.69400
C  -6.98500 22.14200 26.23200
C  -7.94700 20.03100 29.14500
C  -8.63000 20.10000 25.80900
O  -8.90700 20.38500 24.62300
C  -9.08700 18.72400 26.28000
N  -5.32000 24.55800 25.17000
C  -5.97300 23.61500 24.38200
C  -5.69900 23.91000 22.94300
C  -4.77200 25.17400 22.83700
C  -4.60800 25.45200 24.34000
C  -5.24900 22.67900 22.10000
C  -5.37300 26.38000 22.08300
C  -4.94400 26.45700 20.60500
N  -3.89500 26.01800 27.18100
C  -3.38100 26.77900 26.13200
C  -2.50100 27.90900 26.61000
C  -2.55000 27.71500 28.00000
C  -3.42600 26.62300 28.30100
C  -1.78700 28.93500 25.74900
C  -2.03300 28.27400 29.24000
O  -1.24900 29.20200 29.44300
C  -2.78100 27.47500 30.39000
C  -1.66500 26.98400 31.33300
O  -0.58100 26.54900 31.00300
O  -2.03300 27.27600 32.64500
C  -1.03300 26.85400 33.62600
H  -6.79600 21.92300 30.55500
H  -7.24600 22.04300 24.05300
H  -3.20100 26.86300 23.98700
H  -3.94200 25.15700 32.46000
H  -6.32400 23.59600 32.27400
H  -3.58700 22.93500 32.90600
H  -4.09000 21.81900 31.73300
H  -5.03500 22.07100 33.19200
H  -6.83700 25.75100 32.17600
H  -5.77900 27.08800 31.73500
H  -6.23400 27.46900 33.98000
H  -4.65000 26.93300 33.94100
H  -7.01100 19.85400 29.67500
H  -8.29300 19.05500 28.80500
H  -8.67800 20.40400 29.86200
H  -8.25700 18.16700 26.71500
H  -9.45700 18.28200 25.35500
H  -9.85300 18.80700 27.05100
H  -6.66900 24.06000 22.46800
H  -3.77200 24.98400 22.44500
H  -5.16100 21.82900 22.77700
H  -4.26700 22.88500 21.67500
H  -6.05700 22.54500 21.38000
H  -5.16500 27.26900 22.67800
H  -6.45400 26.27700 21.98600
H  -4.65800 25.46600 20.25200
H  -4.05700 27.05100 20.38400
H  -5.81000 26.77000 20.02100
H  -2.54700 29.61400 25.36200
H  -1.26100 28.53800 24.88100
H  -1.10500 29.45700 26.42000
H  -3.50500 28.04800 30.96900
H  -1.16000 25.83800 34.00100
H  -1.12400 27.59400 34.42100
H  0.01800 26.83100 33.33800


