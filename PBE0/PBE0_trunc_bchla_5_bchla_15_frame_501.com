%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg 47.33900 34.91800 28.07800
C  45.68800 33.19800 30.68300
C  48.10100 37.44900 30.24400
C  48.24300 36.72400 25.47200
C  46.07200 32.47700 25.83700
N  47.02000 35.19800 30.25900
C  46.36000 34.38500 31.06000
C  46.67100 34.77400 32.49400
C  47.03300 36.30500 32.38100
C  47.41200 36.35800 30.85400
C  45.84200 37.21500 32.86800
C  47.88800 34.00400 33.06800
C  47.98300 33.94800 34.65600
H  46.95600 34.64300 35.55400
N  48.24300 36.70300 27.90700
C  48.56900 37.63000 28.93000
C  49.30300 38.76000 28.31400
C  49.26700 38.57700 26.88500
C  48.57000 37.29400 26.64900
C  49.93800 39.85800 29.09200
C  49.83100 39.42600 25.68000
O  49.78800 39.06500 24.49600
C  50.32700 40.81500 25.96000
N  47.07800 34.69300 25.91700
C  47.58500 35.54400 25.06500
C  47.43000 35.17200 23.61700
C  46.71000 33.78100 23.71200
C  46.68600 33.56800 25.25900
C  46.76800 36.29200 22.80900
C  47.43400 32.63900 23.01000
C  46.53900 31.60900 22.23200
N  46.16500 33.09500 28.16200
C  45.77400 32.19200 27.15400
C  45.18200 31.05300 27.76700
C  44.95200 31.44600 29.12400
C  45.59700 32.64100 29.31900
C  44.92100 29.73400 27.09800
C  44.48300 31.08700 30.48400
O  43.89600 30.09800 30.87100
C  44.73800 32.34900 31.43800
C  44.99200 31.90700 32.84200
O  45.89000 31.21500 33.22200
O  43.96300 32.32800 33.70500
C  43.67300 31.45600 34.86300
H  48.20800 38.37100 30.82000
H  48.63100 37.17100 24.55500
H  45.76300 31.70600 25.12800
H  45.78600 34.43600 33.03300
H  47.91000 36.55500 32.97900
H  44.99600 36.59900 33.17100
H  45.50600 37.70100 31.95200
H  46.07700 37.86100 33.71400
H  48.75200 34.55000 32.68700
H  47.86600 33.00400 32.63500
H  48.97800 34.33500 34.88100
H  48.01300 32.88400 34.89200
H  50.02200 39.76300 30.17400
H  49.40400 40.80200 28.98700
H  50.95300 39.89900 28.69600
H  51.37600 40.86200 26.25200
H  49.62900 41.09200 26.75100
H  50.11900 41.39300 25.06000
H  48.46800 35.04700 23.30900
H  45.69100 33.78300 23.32700
H  46.13600 36.95200 23.40400
H  46.11300 35.78300 22.10300
H  47.49100 36.91500 22.28200
H  48.02000 32.05200 23.71700
H  48.21400 33.13900 22.43600
H  45.47900 31.78000 22.41900
H  46.81500 30.64200 22.65200
H  46.65400 31.52300 21.15100
H  44.41500 29.99000 26.16700
H  44.28000 29.05700 27.66300
H  45.86600 29.26100 26.82900
H  43.79600 32.88200 31.56300
H  43.90900 32.01600 35.76800
H  44.22800 30.51800 34.82500
H  42.60300 31.25600 34.80700


