%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 44.01800 2.91800 47.66100
C  42.06000 5.71500 47.01900
C  41.00600 1.01600 47.31200
C  45.79900 0.15700 47.76200
C  46.87900 4.91700 47.24100
N  41.83200 3.33000 47.09800
C  41.29400 4.59200 46.98100
C  39.79500 4.50000 46.68300
C  39.48700 3.02600 47.06600
C  40.84700 2.38200 47.12900
C  38.69200 2.74400 48.40100
C  39.44800 4.83900 45.22400
C  38.00300 4.93400 44.91400
H  37.59500 5.76800 43.64000
N  43.49300 0.83400 47.54100
C  42.21400 0.26700 47.42800
C  42.35900 -1.16200 47.43400
C  43.77600 -1.43800 47.42700
C  44.43200 -0.09700 47.60800
C  41.15600 -2.13400 47.28600
C  44.34200 -2.85900 47.34100
O  43.62700 -3.84400 47.42200
C  45.80100 -3.04100 47.16300
N  46.06200 2.62900 47.74100
C  46.53800 1.38200 47.82400
C  48.11500 1.40300 47.69600
C  48.51000 2.91300 47.47000
C  47.02000 3.56100 47.42600
C  48.93400 0.68700 48.80200
C  49.24600 3.07700 46.12300
C  50.25900 4.22700 46.06800
N  44.47400 4.92000 47.22100
C  45.65200 5.61300 47.20700
C  45.34700 7.04400 47.09000
C  43.98000 7.11700 47.10900
C  43.47900 5.83500 47.12000
C  46.36600 8.18800 47.10600
C  42.87400 8.02300 47.02800
O  42.81600 9.22700 46.97900
C  41.57700 7.14800 46.98100
C  40.79100 7.44500 48.21000
O  41.19100 7.46200 49.34100
O  39.43700 7.57100 47.85400
C  38.52600 7.78100 48.94300
H  40.02600 0.53600 47.28500
H  46.51600 -0.66500 47.82200
H  47.74600 5.57800 47.17700
H  39.30800 5.14400 47.41400
H  38.96600 2.43700 46.31200
H  37.88000 2.08800 48.08800
H  38.53600 3.57700 49.08700
H  39.33000 2.09900 49.00500
H  39.75400 4.10000 44.48300
H  39.93900 5.76800 44.93400
H  37.52300 5.35800 45.79600
H  37.59600 3.92700 44.82400
H  41.34100 -2.41600 46.24900
H  40.16900 -1.67200 47.25600
H  41.01800 -3.02000 47.90600
H  46.05900 -2.38200 46.33400
H  45.98200 -4.06600 46.83900
H  46.21200 -2.73700 48.12600
H  48.33800 0.83400 46.79300
H  49.07900 3.37000 48.27900
H  49.32800 -0.23900 48.38200
H  48.30000 0.56800 49.68000
H  49.82500 1.27000 49.03500
H  48.45300 3.22300 45.39000
H  49.73200 2.16800 45.76800
H  51.18600 3.72000 45.79800
H  50.32400 4.81900 46.98100
H  50.02100 4.90100 45.24500
H  46.22000 8.97700 46.36900
H  47.30400 7.65100 46.96300
H  46.44100 8.50600 48.14600
H  41.03900 7.44900 46.08200
H  38.18200 8.81200 48.86000
H  38.86300 7.60600 49.96500
H  37.69400 7.09100 48.80500
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


