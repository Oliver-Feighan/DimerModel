%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 46.26900 25.13700 28.87000
C  47.26400 27.26100 31.45200
C  45.06800 22.98000 31.10200
C  46.59300 22.69600 26.49900
C  47.58600 27.40300 26.57500
N  46.02800 25.20400 31.05400
C  46.48600 26.15500 31.85900
C  45.98900 25.86800 33.32200
C  45.39300 24.44600 33.19000
C  45.51200 24.15500 31.69700
C  46.24000 23.30800 33.94900
C  44.92100 26.90700 33.89400
C  45.18000 27.25800 35.31900
H  44.05100 27.15600 36.33300
N  45.85600 23.12600 28.82000
C  45.33200 22.40200 29.85800
C  44.87200 21.10600 29.34000
C  45.32500 21.03600 27.99300
C  46.00600 22.31800 27.74700
C  44.15900 20.03600 30.22700
C  45.21800 19.84800 26.99000
O  45.63600 19.88000 25.83800
C  44.54500 18.61000 27.44900
N  46.98900 25.03900 26.90500
C  47.10000 23.90200 26.09900
C  47.90400 24.25000 24.81900
C  47.61600 25.80400 24.63200
C  47.36900 26.07600 26.14300
C  49.38900 23.83600 24.90100
C  46.42100 26.18700 23.80600
C  45.03400 26.12100 24.37300
N  47.05200 27.04800 28.89500
C  47.58500 27.84200 27.90600
C  48.07800 29.11700 28.43200
C  48.07000 28.90600 29.88400
C  47.43900 27.61600 30.07900
C  48.67900 30.18600 27.67300
C  48.43700 29.46600 31.17900
O  49.02500 30.43400 31.48700
C  47.84100 28.42700 32.26500
C  48.81200 28.03800 33.35100
O  49.55700 27.06300 33.31900
O  48.64800 28.91100 34.42600
C  49.48500 28.68900 35.61900
H  44.47200 22.24900 31.65300
H  46.68200 21.86300 25.79800
H  47.89200 28.10500 25.79700
H  46.81500 25.79600 34.03000
H  44.37900 24.31000 33.56600
H  45.81000 23.21000 34.94600
H  47.29000 23.54400 34.12000
H  46.10500 22.37300 33.40400
H  43.90600 26.51300 33.84000
H  44.83700 27.84500 33.34500
H  45.39600 28.32500 35.27200
H  46.03100 26.79100 35.81400
H  44.95600 19.41500 30.63500
H  43.40600 19.48400 29.66400
H  43.77200 20.42600 31.16900
H  44.70200 18.21100 28.45100
H  44.79100 17.78900 26.77600
H  43.48800 18.66000 27.18600
H  47.38400 23.69600 24.03800
H  48.47400 26.33800 24.22400
H  50.12500 24.61500 24.70200
H  49.61300 22.90600 24.37800
H  49.63200 23.53200 25.91900
H  46.32100 25.62700 22.87700
H  46.56100 27.21900 23.48400
H  44.87400 27.14600 24.70700
H  44.90600 25.48200 25.24700
H  44.33600 25.84200 23.58300
H  48.36400 31.10800 28.16200
H  48.25000 30.14800 26.67200
H  49.76500 30.10200 27.69900
H  47.00200 28.92500 32.75100
H  50.51900 28.64400 35.27600
H  49.25400 27.73000 36.08200
H  49.44300 29.39800 36.44600
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


