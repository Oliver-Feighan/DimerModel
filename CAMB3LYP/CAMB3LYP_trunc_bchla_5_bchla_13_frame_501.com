%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 24.36900 -6.98700 46.60800
C  26.68400 -4.76900 45.35200
C  21.90100 -5.12400 45.15800
C  22.46500 -9.68400 46.72500
C  27.03400 -8.80400 47.99500
N  24.41600 -5.40800 45.00500
C  25.43000 -4.51800 44.85100
C  25.02900 -3.31400 44.05100
C  23.50400 -3.26400 44.39000
C  23.25600 -4.65400 44.95700
C  23.12700 -2.08400 45.33100
C  25.27600 -3.52100 42.47300
C  25.72300 -2.21200 41.67100
H  24.68300 -1.93800 40.61300
N  22.38400 -7.27700 46.19700
C  21.52700 -6.46900 45.53500
C  20.22900 -7.11300 45.37800
C  20.37200 -8.43000 45.86100
C  21.79800 -8.50600 46.28600
C  19.04300 -6.37400 44.71200
C  19.31200 -9.48200 45.68900
O  18.25700 -9.18500 45.19200
C  19.48300 -10.96700 46.06200
N  24.62900 -8.97800 47.39100
C  23.73100 -9.85500 47.30000
C  24.17900 -11.20600 47.72000
C  25.67400 -10.94700 48.15800
C  25.81200 -9.52100 47.84300
C  23.22200 -11.82700 48.79300
C  26.70500 -11.85900 47.34500
C  28.13900 -12.05100 47.95500
N  26.42100 -6.83200 46.76500
C  27.35000 -7.59400 47.45000
C  28.63900 -6.96000 47.42800
C  28.41400 -5.78200 46.63300
C  27.10800 -5.81900 46.18400
C  29.90100 -7.42400 48.08100
C  29.04200 -4.60800 46.13300
O  30.18200 -4.13000 46.19600
C  28.00800 -4.00400 45.12600
C  27.97900 -2.53300 45.37600
O  28.65400 -1.69700 44.75900
O  27.14200 -2.25800 46.42200
C  27.16400 -0.79200 46.68700
H  21.02900 -4.53900 44.85700
H  21.83100 -10.57200 46.69600
H  27.87200 -9.43200 48.30700
H  25.56100 -2.49300 44.53000
H  22.94300 -3.07800 43.47400
H  22.46400 -1.40000 44.80200
H  24.04300 -1.52500 45.52500
H  22.72200 -2.51300 46.24700
H  24.47600 -4.01300 41.92000
H  26.07600 -4.24200 42.30800
H  26.63900 -2.51400 41.16300
H  25.99300 -1.30800 42.21700
H  18.63800 -6.86100 43.82500
H  19.24900 -5.33100 44.46900
H  18.26300 -6.29000 45.46800
H  19.62500 -11.12900 47.13000
H  20.27800 -11.46900 45.51000
H  18.60200 -11.45500 45.64400
H  24.15900 -11.77400 46.79000
H  25.92100 -11.00100 49.21800
H  22.84600 -12.76100 48.37600
H  22.40400 -11.15300 49.04800
H  23.73300 -12.22400 49.67000
H  26.84400 -11.46000 46.34000
H  26.21500 -12.82700 47.25100
H  28.51200 -13.04100 48.21600
H  28.22300 -11.52500 48.90500
H  28.89100 -11.63100 47.28700
H  29.81500 -8.26400 48.77000
H  30.19100 -6.49700 48.57400
H  30.52400 -7.56700 47.19800
H  28.23100 -4.20700 44.07800
H  28.13300 -0.51500 47.10100
H  26.50700 -0.50200 47.50700
H  26.90400 -0.28600 45.75700
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


