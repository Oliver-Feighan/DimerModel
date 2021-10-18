%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 8.37300 3.02000 27.58300
C  9.66800 1.59900 30.62700
C  6.70500 5.26800 29.53700
C  6.93600 4.14100 24.76700
C  9.73100 0.40000 25.89700
N  8.26000 3.41400 29.92600
C  8.92300 2.70600 30.88600
C  8.72800 3.42100 32.16900
C  7.67000 4.51900 31.87000
C  7.43900 4.42400 30.33700
C  6.33800 4.27100 32.64700
C  9.97700 3.96400 32.81700
C  9.90300 4.49900 34.23800
H  10.91500 4.00800 35.18600
N  7.11500 4.54900 27.20100
C  6.43400 5.32900 28.17600
C  5.69100 6.36700 27.43700
C  5.70300 5.96800 26.10300
C  6.57300 4.78400 25.98300
C  4.94900 7.46900 28.16400
C  5.04400 6.62700 24.87200
O  4.92600 6.07300 23.76300
C  4.55800 8.00900 25.03700
N  8.26900 2.33700 25.58000
C  7.65300 2.93000 24.57800
C  7.79400 2.17100 23.23000
C  8.87800 1.15300 23.59100
C  9.00200 1.33000 25.08700
C  6.44500 1.50000 22.80300
C  10.20800 1.29100 22.88500
C  10.86600 -0.01000 22.39200
N  9.53800 1.41400 28.13000
C  10.02500 0.42800 27.30000
C  10.88100 -0.54000 28.05000
C  10.69600 -0.06300 29.33200
C  9.85900 1.05900 29.36100
C  11.78700 -1.59900 27.50900
C  11.07900 -0.33300 30.71000
O  11.83900 -1.20400 31.11800
C  10.43500 0.73200 31.62700
C  9.77600 0.01400 32.68800
O  8.68700 -0.60500 32.61700
O  10.41500 0.23000 33.87100
C  9.76200 -0.17400 35.09700
H  6.23400 6.04600 30.14100
H  6.34700 4.39700 23.88400
H  9.91800 -0.50400 25.31400
H  8.28800 2.65400 32.80500
H  8.01700 5.51600 32.14200
H  5.84500 5.03100 33.25400
H  6.57700 3.45100 33.32400
H  5.53900 3.81600 32.06100
H  10.35300 4.83500 32.28000
H  10.84000 3.32300 32.63900
H  8.94300 4.25600 34.69100
H  10.21000 5.54200 34.16500
H  3.93400 7.63800 27.80500
H  5.45900 8.43200 28.14900
H  4.76300 7.27700 29.22100
H  5.36300 8.60800 25.46300
H  3.75000 7.88400 25.75700
H  4.37900 8.39200 24.03300
H  8.13900 2.90800 22.50400
H  8.61200 0.10400 23.46400
H  6.50600 0.41300 22.84700
H  6.19800 1.73200 21.76700
H  5.70700 1.85800 23.52100
H  10.93600 1.82600 23.49400
H  10.03500 2.08200 22.15600
H  11.47800 0.25700 21.53100
H  10.15100 -0.81800 22.23400
H  11.45900 -0.28500 23.26400
H  12.83900 -1.31500 27.47600
H  11.33000 -1.91300 26.57100
H  11.63300 -2.52100 28.07100
H  11.28200 1.20800 32.12100
H  8.73800 -0.52700 34.97200
H  9.83700 0.61100 35.84900
H  10.27900 -1.00300 35.58100
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


