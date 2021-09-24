%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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


