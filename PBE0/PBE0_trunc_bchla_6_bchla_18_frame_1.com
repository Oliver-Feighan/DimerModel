%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 17.05400 -2.28800 27.75000
C  16.54800 -0.35400 30.72600
C  18.85200 -4.52900 29.65700
C  18.01300 -3.67000 24.94900
C  15.39200 0.22200 26.03300
N  17.88600 -2.27700 29.96200
C  17.41200 -1.43700 30.93300
C  17.78200 -2.00100 32.33100
C  18.58500 -3.32700 31.99600
C  18.48400 -3.42400 30.43000
C  20.10200 -3.41100 32.45300
C  16.57800 -2.28600 33.26800
C  16.68400 -1.85000 34.68600
H  18.05300 -1.52300 35.39600
N  18.30500 -3.91800 27.31000
C  18.85900 -4.81900 28.21900
C  19.64100 -5.82500 27.47500
C  19.29300 -5.63200 26.08300
C  18.51600 -4.37100 26.05800
C  20.53100 -6.87100 28.07500
C  19.63900 -6.46300 24.82500
O  19.04500 -6.38400 23.72900
C  20.63000 -7.54600 24.97700
N  16.79100 -1.75400 25.84400
C  17.37300 -2.49900 24.84200
C  16.84500 -2.07000 23.44500
C  15.70300 -1.04100 23.86200
C  15.95300 -0.84300 25.34700
C  17.88000 -1.47900 22.40300
C  14.26500 -1.46100 23.60200
C  13.50200 -0.61900 22.55800
N  15.98100 -0.53800 28.25500
C  15.38600 0.39400 27.43800
C  14.87800 1.51400 28.24100
C  15.27300 1.20700 29.53100
C  15.93800 -0.03900 29.49700
C  14.17200 2.75100 27.79300
C  15.18900 1.69900 30.85900
O  14.57400 2.68200 31.36200
C  15.98100 0.75800 31.66400
C  17.01200 1.50300 32.46800
O  17.97500 2.03100 32.02700
O  16.72000 1.42300 33.76600
C  17.36800 2.29700 34.68400
H  19.44500 -5.28100 30.18200
H  18.21500 -4.04600 23.94300
H  14.94000 1.04300 25.47400
H  18.42200 -1.24900 32.79400
H  18.07000 -4.18900 32.41900
H  20.43900 -2.60400 33.10300
H  20.70600 -3.31100 31.55000
H  20.31100 -4.38200 32.90300
H  16.34900 -3.34900 33.35200
H  15.66000 -1.80400 32.93300
H  16.30900 -2.63700 35.34100
H  16.04900 -0.97200 34.80100
H  21.52400 -6.79800 27.63200
H  20.09600 -7.83300 27.80400
H  20.47400 -6.87000 29.16400
H  21.57900 -7.08500 25.25100
H  20.75500 -7.96300 23.97800
H  20.26200 -8.35300 25.61000
H  16.37000 -2.96000 23.03100
H  15.88800 -0.12500 23.30100
H  18.88000 -1.46500 22.83500
H  17.57100 -0.51800 21.99000
H  17.86400 -2.13300 21.53100
H  13.63500 -1.61600 24.47800
H  14.30700 -2.45900 23.16700
H  12.47900 -0.33500 22.80400
H  13.29800 -1.28100 21.71600
H  13.99200 0.30400 22.25100
H  13.75500 2.65600 26.79000
H  14.88900 3.56000 27.65400
H  13.49700 3.03000 28.60200
H  15.27200 0.27100 32.33300
H  16.91900 3.28000 34.53500
H  18.41700 2.43200 34.42000
H  17.17400 2.00000 35.71500
Mg 35.49600 50.25400 25.27000
C  35.17700 48.23300 28.13300
C  33.95100 52.80200 27.06200
C  35.52400 52.00900 22.45100
C  36.25500 47.31600 23.35100
N  34.65000 50.47900 27.32000
C  34.75900 49.53200 28.34800
C  34.55800 50.22300 29.68400
C  33.90100 51.58900 29.33600
C  34.18900 51.63300 27.80500
C  32.41200 51.77300 29.76700
C  35.85100 50.33500 30.56400
C  35.57700 50.30300 32.13200
H  34.12200 50.22500 32.60200
N  34.96300 52.12400 24.87400
C  34.41800 53.04500 25.75200
C  34.29200 54.33000 25.12700
C  34.61900 54.10800 23.75100
C  35.08400 52.73500 23.66700
C  33.72300 55.59000 25.84800
C  34.52600 55.05800 22.60200
O  34.69900 54.67600 21.42000
C  34.14200 56.56700 22.60300
N  35.57000 49.71100 23.15300
C  35.66800 50.65100 22.19400
C  36.16600 50.02300 20.86800
C  36.37900 48.52900 21.13700
C  36.08600 48.52400 22.70200
C  35.30200 50.41100 19.70700
C  37.75400 47.93500 20.76300
C  38.97000 48.58500 21.46800
N  35.63500 48.14500 25.59300
C  36.09300 47.16800 24.77600
C  36.32500 45.99400 25.61600
C  36.01700 46.32200 26.96600
C  35.54300 47.65500 26.90400
C  36.93800 44.67500 25.15300
C  35.96700 45.83900 28.30900
O  36.34100 44.80600 28.85800
C  35.32300 47.11200 29.16700
C  34.10500 46.63100 29.93300
O  33.00700 46.49100 29.43100
O  34.39500 46.48000 31.27000
C  33.50200 45.67800 32.10000
H  33.57900 53.56200 27.75200
H  35.60400 52.66100 21.57900
H  36.61900 46.49900 22.72500
H  33.85700 49.64500 30.28600
H  34.49900 52.36000 29.82100
H  31.92700 52.02300 28.82400
H  32.30200 52.65000 30.40600
H  31.87100 50.93900 30.21400
H  36.34800 51.29400 30.42100
H  36.54700 49.55500 30.25800
H  35.89400 51.26500 32.53500
H  36.14900 49.47100 32.54200
H  33.11100 56.13200 25.12600
H  34.50500 56.27400 26.17600
H  33.12800 55.50300 26.75700
H  34.79000 57.12800 23.27700
H  33.09000 56.73000 22.83700
H  34.45500 57.11000 21.71200
H  37.07600 50.56600 20.61600
H  35.49800 48.04500 20.71400
H  35.74600 51.19400 19.09200
H  34.31200 50.69900 20.06100
H  35.07200 49.59000 19.02700
H  37.83100 48.14500 19.69600
H  37.75500 46.85700 20.92700
H  39.44400 47.77500 22.02200
H  38.57400 49.35600 22.12900
H  39.64500 49.11900 20.79900
H  37.67500 44.28600 25.85500
H  37.52900 44.80600 24.24700
H  36.12400 43.98800 24.92000
H  36.09300 47.35400 29.90000
H  32.48700 45.91500 31.78000
H  33.56400 45.88200 33.16900
H  33.81500 44.64500 31.94800


