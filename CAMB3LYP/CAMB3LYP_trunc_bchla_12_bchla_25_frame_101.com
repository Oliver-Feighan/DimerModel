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
Mg -2.59000 34.01700 26.79400
C  -3.62000 32.41400 29.77300
C  -1.15300 36.45800 28.71900
C  -2.20600 35.82500 24.04600
C  -4.51700 31.70900 25.00200
N  -2.55200 34.47200 29.00300
C  -2.95200 33.59300 30.01300
C  -2.37700 34.05700 31.37400
C  -1.60200 35.37500 31.00700
C  -1.71600 35.46600 29.46000
C  -1.91400 36.61400 31.82800
C  -1.50900 33.04200 32.13400
C  -1.98900 32.76400 33.58500
H  -0.95800 32.99700 34.66700
N  -1.70200 35.86000 26.44900
C  -1.19200 36.74900 27.34900
C  -0.55100 37.80700 26.60500
C  -0.98500 37.71000 25.21400
C  -1.68000 36.43200 25.23200
C  0.22600 39.02700 27.32500
C  -1.03600 38.70200 24.13300
O  -1.40600 38.53400 23.00600
C  -0.57600 40.11000 24.49900
N  -3.46600 33.87800 24.82300
C  -3.08300 34.73100 23.85300
C  -3.52500 34.22200 22.47900
C  -3.92700 32.76700 22.76900
C  -3.94500 32.73700 24.28100
C  -4.74000 35.08300 21.90100
C  -2.77000 31.79500 22.20500
C  -1.44300 31.78800 23.01200
N  -3.75200 32.29300 27.19100
C  -4.44100 31.48200 26.37200
C  -5.05500 30.45500 27.17000
C  -4.72100 30.75200 28.47800
C  -3.96900 31.93000 28.47700
C  -5.92100 29.29000 26.73000
C  -4.97000 30.27500 29.84100
O  -5.63000 29.44100 30.39400
C  -4.20300 31.28400 30.73600
C  -5.22200 31.86200 31.64700
O  -6.22400 32.38000 31.22200
O  -4.94100 31.66200 32.98900
C  -6.05600 31.94700 33.95300
H  -0.52600 37.10900 29.33200
H  -2.00600 36.32000 23.09300
H  -4.95700 30.83400 24.51900
H  -3.20300 34.37600 32.01000
H  -0.55500 35.23200 31.27500
H  -2.84600 36.49700 32.37900
H  -1.95200 37.53200 31.24100
H  -1.13500 36.66900 32.58800
H  -0.45800 33.28400 32.29200
H  -1.47400 32.09100 31.60300
H  -2.46500 31.78700 33.66800
H  -2.77300 33.48500 33.81300
H  1.18600 39.15100 26.82400
H  0.47900 38.88500 28.37500
H  -0.26700 39.99700 27.38900
H  0.46900 39.96700 24.77500
H  -1.19700 40.39100 25.34900
H  -0.59300 40.82400 23.67500
H  -2.71600 34.17900 21.74900
H  -4.90200 32.39200 22.45700
H  -4.66500 35.04300 20.81400
H  -4.79700 36.14200 22.15500
H  -5.62600 34.53900 22.22800
H  -2.53100 31.98600 21.15900
H  -3.07000 30.75600 22.34200
H  -0.61000 31.26900 22.53900
H  -1.63900 31.28300 23.95800
H  -1.17900 32.82600 23.21700
H  -5.41600 28.57400 26.08200
H  -6.72000 29.76700 26.16300
H  -6.40800 28.78800 27.56600
H  -3.52700 30.64600 31.30500
H  -5.48500 32.48600 34.70900
H  -6.48400 31.04500 34.39000
H  -6.92400 32.52400 33.63500


