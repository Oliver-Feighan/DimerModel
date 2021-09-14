%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 46.63400 24.73200 28.52400
C  47.16600 26.79400 31.26900
C  45.35000 22.41900 30.76700
C  46.57500 22.51000 25.95000
C  48.10800 27.10900 26.49100
N  46.35100 24.60700 30.77200
C  46.60600 25.59200 31.71400
C  46.32600 25.14000 33.20100
C  46.05100 23.60400 32.96700
C  45.94200 23.50100 31.40100
C  47.03100 22.54000 33.47000
C  45.17300 25.95500 33.81500
C  45.53600 26.98000 34.85700
H  44.55500 27.15400 36.03000
N  45.88000 22.81500 28.28900
C  45.34100 22.11500 29.32300
C  44.95500 20.80900 28.75400
C  45.37600 20.75700 27.41600
C  45.95000 22.01600 27.13800
C  44.14400 19.85300 29.53100
C  45.30700 19.63400 26.42300
O  45.62600 19.75900 25.20800
C  44.90700 18.30200 26.75200
N  47.16200 24.82000 26.54200
C  47.02200 23.79200 25.66100
C  47.50000 24.11000 24.29900
C  47.94400 25.61100 24.43200
C  47.78600 25.82800 25.89000
C  48.58100 23.27100 23.55600
C  47.18600 26.64400 23.46000
C  45.80900 27.28400 23.80700
N  47.51200 26.55300 28.75500
C  48.03300 27.42300 27.84100
C  48.30100 28.65400 28.44600
C  48.04200 28.39000 29.84600
C  47.54000 27.10900 29.96100
C  48.74800 29.95300 27.77800
C  47.98400 29.05500 31.11600
O  48.14900 30.19900 31.48800
C  47.61900 27.93100 32.14800
C  48.80500 27.58900 32.95800
O  49.43200 26.50700 32.93800
O  49.23500 28.58100 33.76800
C  50.27600 28.29000 34.81200
H  45.00700 21.67200 31.48500
H  46.61700 21.83300 25.09400
H  48.42100 27.96800 25.89300
H  47.24300 25.06200 33.78600
H  45.06600 23.41800 33.39500
H  47.53600 22.06600 32.62800
H  46.41900 21.72900 33.86200
H  47.76900 22.88500 34.19500
H  44.55000 25.22000 34.32600
H  44.58000 26.37000 33.00000
H  45.61600 27.94300 34.35200
H  46.51600 26.85600 35.31800
H  43.21600 19.75600 28.96900
H  43.90200 20.07100 30.57200
H  44.63900 18.88200 29.54400
H  45.42400 18.07100 27.68300
H  45.01300 17.64200 25.89200
H  43.86100 18.27000 27.05700
H  46.55600 24.02200 23.76200
H  49.01300 25.68800 24.23700
H  48.17600 22.76600 22.67900
H  48.97300 22.47900 24.19500
H  49.46000 23.85000 23.27300
H  47.00900 26.03100 22.57600
H  47.98400 27.36400 23.27800
H  45.45700 27.11300 24.82400
H  45.01200 27.13000 23.07900
H  45.99000 28.35100 23.68400
H  49.68100 30.32800 28.20000
H  47.99600 30.71000 28.00000
H  48.80300 29.90800 26.69000
H  46.94500 28.37300 32.88200
H  50.12900 29.19800 35.39700
H  51.29200 28.18100 34.43400
H  50.01200 27.37900 35.35000
Mg -5.09800 24.63300 26.68400
C  -3.70200 26.36700 29.52200
C  -6.11300 22.14300 28.93500
C  -6.48000 23.09600 24.13300
C  -3.72900 27.06700 24.58100
N  -4.87100 24.33000 29.01500
C  -4.47600 25.23800 29.94600
C  -4.87000 24.81900 31.35300
C  -5.34300 23.26900 31.18300
C  -5.42800 23.22100 29.63000
C  -4.42000 22.22500 31.86000
C  -5.97200 25.74300 31.97300
C  -5.56400 26.39400 33.32000
H  -6.48200 25.91600 34.58700
N  -6.12400 22.92100 26.55900
C  -6.45400 21.95500 27.56900
C  -7.23800 20.80800 27.05300
C  -7.37200 21.15200 25.66800
C  -6.65700 22.41400 25.41400
C  -7.70400 19.61500 27.89000
C  -8.01000 20.37100 24.58300
O  -8.30600 20.80000 23.47300
C  -8.47200 18.93600 24.90400
N  -5.23600 25.12700 24.71500
C  -5.84300 24.30800 23.82600
C  -5.57400 24.69000 22.39700
C  -4.73400 25.98700 22.49800
C  -4.58300 26.13300 24.02900
C  -4.97900 23.58700 21.54700
C  -5.55400 27.25800 21.86300
C  -4.76800 28.10300 20.72000
N  -3.71900 26.23300 26.92800
C  -3.33300 27.18000 25.94900
C  -2.55400 28.22500 26.62300
C  -2.81300 27.99000 27.99900
C  -3.40400 26.73200 28.12600
C  -1.71900 29.35800 26.03400
C  -2.39200 28.46300 29.34800
O  -1.66000 29.39300 29.78700
C  -3.08100 27.46800 30.36000
C  -2.03500 26.94900 31.30200
O  -0.94600 26.57900 30.96500
O  -2.47500 27.07900 32.59000
C  -1.38400 26.91100 33.58300
H  -6.28600 21.33800 29.65200
H  -6.91400 22.63300 23.24400
H  -3.32400 27.85900 23.94800
H  -3.94900 24.82500 31.93600
H  -6.36100 23.10200 31.53600
H  -3.94100 21.58900 31.11500
H  -4.98800 21.61200 32.56100
H  -3.64700 22.62700 32.51500
H  -6.89500 25.19600 32.16700
H  -6.22100 26.57200 31.31000
H  -5.77000 27.46400 33.35800
H  -4.50500 26.32600 33.56700
H  -7.71300 19.74200 28.97300
H  -7.20000 18.67100 27.68300
H  -8.74200 19.45800 27.59500
H  -7.62000 18.31700 25.18500
H  -9.03200 18.48700 24.08300
H  -9.15600 18.98900 25.75100
H  -6.53500 24.86100 21.91000
H  -3.75100 25.85800 22.04600
H  -5.68000 23.16900 20.82400
H  -4.66700 22.69800 22.09500
H  -4.09700 23.98900 21.04800
H  -5.86900 27.99500 22.60100
H  -6.48500 26.95600 21.38400
H  -4.00100 27.42900 20.33800
H  -4.32800 29.01700 21.12100
H  -5.40200 28.28200 19.85200
H  -1.71900 30.24300 26.67000
H  -2.23000 29.70200 25.13400
H  -0.77900 28.87000 25.77600
H  -3.74400 28.08300 30.96800
H  -0.90300 25.95500 33.37700
H  -1.91800 26.79500 34.52600
H  -0.58100 27.64700 33.62200


