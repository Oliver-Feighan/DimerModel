%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 40.74100 8.40300 29.69900
C  42.66000 10.02500 32.01600
C  38.73200 7.21500 32.07700
C  39.17900 6.75200 27.36300
C  43.28300 9.22800 27.34100
N  40.79300 8.45600 31.86800
C  41.50500 9.39200 32.54400
C  40.94800 9.63500 33.92000
C  39.85300 8.49700 34.00200
C  39.82700 7.94800 32.59200
C  40.10000 7.44100 35.10600
C  40.44100 11.18100 34.18500
C  40.87000 11.88200 35.48300
H  41.23900 10.96600 36.69200
N  39.13800 7.34400 29.69200
C  38.39900 6.91700 30.74900
C  37.34700 6.00300 30.37600
C  37.33700 6.02300 29.01400
C  38.65500 6.62800 28.64700
C  36.61000 5.29800 31.50000
C  36.33400 5.36600 28.04500
O  36.44200 5.38700 26.84100
C  35.20700 4.65100 28.73000
N  41.30000 7.91900 27.67200
C  40.36700 7.31400 26.91500
C  40.80200 7.26400 25.44700
C  41.95100 8.32700 25.38700
C  42.23200 8.49800 26.89400
C  41.21600 5.86100 24.82900
C  41.58000 9.72400 24.76100
C  40.48400 10.59900 25.31000
N  42.57200 9.48100 29.57000
C  43.55300 9.66900 28.59100
C  44.67900 10.46300 29.09900
C  44.36400 10.63700 30.45500
C  43.08700 10.02200 30.64800
C  45.85600 11.01300 28.33500
C  44.78800 11.17000 31.66800
O  45.80500 11.75000 32.01500
C  43.61600 11.01200 32.67100
C  44.29700 10.49000 33.87800
O  44.63700 9.31800 34.05700
O  44.61300 11.49900 34.74100
C  45.54500 11.24000 35.84300
H  37.99000 7.00700 32.85100
H  38.68100 6.27000 26.52000
H  43.95700 9.44200 26.50800
H  41.68000 9.37000 34.68300
H  38.87700 8.95900 34.14900
H  40.48500 6.53200 34.64400
H  39.14600 7.20600 35.57800
H  40.83000 7.87400 35.79000
H  39.35300 11.12200 34.15300
H  40.80700 11.84200 33.40000
H  40.04300 12.46600 35.88600
H  41.65400 12.60800 35.26800
H  37.22900 4.66200 32.13300
H  35.89000 4.51000 31.28000
H  36.12800 6.14400 31.99000
H  34.67500 5.32500 29.40100
H  35.72300 3.80100 29.17600
H  34.37800 4.27900 28.12900
H  39.88500 7.56300 24.93800
H  42.83200 7.86300 24.94300
H  41.11200 5.08600 25.58800
H  42.28100 5.83900 24.60100
H  40.70000 5.56000 23.91600
H  41.36300 9.40900 23.74000
H  42.43100 10.40400 24.73900
H  39.92000 10.04700 26.06300
H  39.77700 10.86400 24.52400
H  40.79900 11.50000 25.83600
H  46.78200 10.48700 28.56300
H  45.97200 12.07500 28.55200
H  45.60900 10.77100 27.30100
H  43.07400 11.94700 32.81400
H  46.05900 12.14500 36.16900
H  46.18500 10.38200 35.64300
H  45.04400 10.86200 36.73400
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


