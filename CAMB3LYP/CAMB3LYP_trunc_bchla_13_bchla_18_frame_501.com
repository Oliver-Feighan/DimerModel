%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg 35.28300 49.28100 25.05100
C  35.19500 47.66600 28.21300
C  33.70700 52.00200 26.67100
C  34.96700 50.74000 22.15000
C  36.56200 46.36300 23.65600
N  34.51200 49.78400 27.24100
C  34.71400 48.96400 28.32400
C  34.31200 49.60400 29.64200
C  33.71300 50.96400 29.08000
C  34.02200 50.95100 27.57500
C  32.25500 51.13300 29.49700
C  35.54400 49.80300 30.55200
C  35.36900 49.62200 32.10800
H  33.94300 49.72800 32.69800
N  34.52100 51.15700 24.48200
C  33.85000 52.07800 25.28700
C  33.47900 53.23000 24.45000
C  33.93900 52.90600 23.13200
C  34.54200 51.56500 23.25500
C  32.73100 54.35200 24.94000
C  33.82500 53.73400 21.85600
O  34.18700 53.27500 20.76100
C  33.44700 55.16200 21.91300
N  35.61000 48.60100 23.15100
C  35.50000 49.43800 22.11600
C  35.83700 48.71600 20.77100
C  36.52000 47.40900 21.28500
C  36.22400 47.41000 22.75200
C  34.60600 48.52300 19.81400
C  38.05900 47.61300 21.01100
C  38.99100 48.19700 22.12800
N  35.88400 47.43300 25.74900
C  36.37900 46.32800 25.03300
C  36.56500 45.20400 25.98400
C  36.06300 45.68400 27.18800
C  35.69100 47.00600 27.06500
C  36.87500 43.87200 25.49400
C  35.71300 45.23100 28.53200
O  35.83100 44.10500 29.00200
C  35.09700 46.50900 29.27300
C  33.71900 46.24000 29.75100
O  32.68000 46.50000 29.15100
O  33.78400 45.64900 31.02700
C  32.59200 45.22000 31.72300
H  33.45300 52.96300 27.12400
H  34.90900 51.17800 21.15100
H  37.00400 45.52600 23.11200
H  33.56100 49.02400 30.17900
H  34.21800 51.75900 29.62800
H  32.03300 52.05800 30.02800
H  31.83200 50.33100 30.10200
H  31.68100 51.22500 28.57400
H  35.98600 50.78700 30.39600
H  36.22100 49.06000 30.12900
H  35.89600 50.41200 32.64200
H  35.97000 48.74800 32.35900
H  33.30100 55.24000 24.66600
H  32.41200 54.18900 25.96900
H  31.80500 54.22200 24.37900
H  32.39300 55.32500 22.14000
H  33.61200 55.61800 20.93700
H  33.97000 55.60700 22.75900
H  36.50600 49.47100 20.35800
H  36.19100 46.53200 20.72700
H  34.58000 49.32300 19.07400
H  33.72600 48.44000 20.45100
H  34.76400 47.56600 19.31700
H  38.19100 48.07200 20.03100
H  38.44400 46.59300 21.03000
H  39.60900 48.91300 21.58700
H  39.65000 47.44800 22.56900
H  38.49600 48.74500 22.92900
H  36.19400 43.63700 24.67700
H  36.66800 43.18600 26.31600
H  37.94500 43.75700 25.32300
H  35.75600 46.63900 30.13200
H  32.37400 44.16000 31.59500
H  31.77700 45.93200 31.59200
H  32.87100 45.44400 32.75300


