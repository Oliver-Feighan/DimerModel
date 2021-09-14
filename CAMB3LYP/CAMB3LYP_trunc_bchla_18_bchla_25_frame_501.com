%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
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
Mg -2.33100 34.27800 27.06700
C  -3.29300 32.73500 29.99200
C  -0.77900 36.82100 29.04000
C  -1.68900 35.99300 24.34500
C  -4.00500 31.86100 25.19100
N  -2.10600 34.64700 29.26400
C  -2.54800 33.87100 30.25800
C  -2.04400 34.45500 31.61800
C  -1.35700 35.84600 31.23800
C  -1.46000 35.80800 29.70100
C  -2.07400 37.11000 31.67200
C  -1.16200 33.50400 32.44100
C  -1.65900 33.38600 33.87200
H  -0.69900 33.03500 35.01800
N  -1.29600 36.10900 26.75300
C  -0.70200 36.96400 27.63600
C  -0.06500 38.02800 26.99800
C  -0.40300 37.90700 25.60500
C  -1.15900 36.61100 25.50300
C  0.74500 39.09500 27.72000
C  -0.12100 38.87200 24.46500
O  -0.65300 38.65500 23.35900
C  0.65800 40.10800 24.66200
N  -2.76200 33.93100 25.11000
C  -2.33900 34.76200 24.09300
C  -2.78500 34.31300 22.65200
C  -3.27700 32.84000 22.98500
C  -3.37800 32.87200 24.52100
C  -3.86900 35.12000 21.83400
C  -2.50400 31.61000 22.32400
C  -1.18700 31.25300 22.96100
N  -3.52800 32.72900 27.43000
C  -4.16000 31.81100 26.61300
C  -4.85300 30.83200 27.34600
C  -4.49800 31.15000 28.66200
C  -3.68000 32.31700 28.68800
C  -5.57400 29.72900 26.69600
C  -4.80700 30.73800 29.97700
O  -5.59800 29.85900 30.32500
C  -4.00900 31.74000 30.96400
C  -4.95800 32.38400 31.85500
O  -5.73200 33.27000 31.57300
O  -4.95600 31.78800 33.04100
C  -5.91200 32.19500 34.02600
H  -0.35300 37.61900 29.65100
H  -1.55300 36.44800 23.36200
H  -4.51800 31.12700 24.56600
H  -2.99900 34.65700 32.10200
H  -0.41300 35.85600 31.78400
H  -2.94800 36.94200 32.30200
H  -2.45800 37.64200 30.80200
H  -1.35600 37.81800 32.08700
H  -0.10200 33.75800 32.45800
H  -1.18900 32.52900 31.95500
H  -2.34600 32.53900 33.87800
H  -2.14000 34.33700 34.10000
H  0.94800 38.81600 28.75400
H  0.13200 39.99100 27.82000
H  1.62000 39.38400 27.13900
H  0.44000 40.83300 23.87700
H  1.70900 39.83900 24.55200
H  0.43200 40.54800 25.63300
H  -1.93500 34.12600 21.99500
H  -4.26000 32.69700 22.53800
H  -4.17100 35.97400 22.43900
H  -4.79700 34.54700 21.82100
H  -3.62200 35.43900 20.82100
H  -2.25200 31.92600 21.31200
H  -3.11500 30.71700 22.19800
H  -1.06000 30.17800 22.83800
H  -0.96900 31.52500 23.99300
H  -0.42200 31.79500 22.40500
H  -5.20900 29.35700 25.73900
H  -6.66000 29.81200 26.67700
H  -5.43800 28.93700 27.43200
H  -3.30800 31.14200 31.54600
H  -5.54400 33.01100 34.64900
H  -6.05600 31.34600 34.69400
H  -6.83100 32.52000 33.53800


