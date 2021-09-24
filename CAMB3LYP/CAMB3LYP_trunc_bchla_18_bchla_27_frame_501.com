%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg -5.09800 24.63300 26.68400
C  -3.70200 26.36700 29.52200
C  -6.11300 22.14300 28.93500
C  -6.48000 23.09600 24.13300
C  -3.72900 27.06700 24.58100
N  -4.87100 24.33000 29.01500
C  -4.47600 25.23800 29.94600
C  -4.87000 24.81900 31.35300
C  -5.34300 23.26900 31.18300
C  -5.42800 23.22100 29.63000
C  -4.42000 22.22500 31.86000
C  -5.97200 25.74300 31.97300
C  -5.56400 26.39400 33.32000
H  -6.48200 25.91600 34.58700
N  -6.12400 22.92100 26.55900
C  -6.45400 21.95500 27.56900
C  -7.23800 20.80800 27.05300
C  -7.37200 21.15200 25.66800
C  -6.65700 22.41400 25.41400
C  -7.70400 19.61500 27.89000
C  -8.01000 20.37100 24.58300
O  -8.30600 20.80000 23.47300
C  -8.47200 18.93600 24.90400
N  -5.23600 25.12700 24.71500
C  -5.84300 24.30800 23.82600
C  -5.57400 24.69000 22.39700
C  -4.73400 25.98700 22.49800
C  -4.58300 26.13300 24.02900
C  -4.97900 23.58700 21.54700
C  -5.55400 27.25800 21.86300
C  -4.76800 28.10300 20.72000
N  -3.71900 26.23300 26.92800
C  -3.33300 27.18000 25.94900
C  -2.55400 28.22500 26.62300
C  -2.81300 27.99000 27.99900
C  -3.40400 26.73200 28.12600
C  -1.71900 29.35800 26.03400
C  -2.39200 28.46300 29.34800
O  -1.66000 29.39300 29.78700
C  -3.08100 27.46800 30.36000
C  -2.03500 26.94900 31.30200
O  -0.94600 26.57900 30.96500
O  -2.47500 27.07900 32.59000
C  -1.38400 26.91100 33.58300
H  -6.28600 21.33800 29.65200
H  -6.91400 22.63300 23.24400
H  -3.32400 27.85900 23.94800
H  -3.94900 24.82500 31.93600
H  -6.36100 23.10200 31.53600
H  -3.94100 21.58900 31.11500
H  -4.98800 21.61200 32.56100
H  -3.64700 22.62700 32.51500
H  -6.89500 25.19600 32.16700
H  -6.22100 26.57200 31.31000
H  -5.77000 27.46400 33.35800
H  -4.50500 26.32600 33.56700
H  -7.71300 19.74200 28.97300
H  -7.20000 18.67100 27.68300
H  -8.74200 19.45800 27.59500
H  -7.62000 18.31700 25.18500
H  -9.03200 18.48700 24.08300
H  -9.15600 18.98900 25.75100
H  -6.53500 24.86100 21.91000
H  -3.75100 25.85800 22.04600
H  -5.68000 23.16900 20.82400
H  -4.66700 22.69800 22.09500
H  -4.09700 23.98900 21.04800
H  -5.86900 27.99500 22.60100
H  -6.48500 26.95600 21.38400
H  -4.00100 27.42900 20.33800
H  -4.32800 29.01700 21.12100
H  -5.40200 28.28200 19.85200
H  -1.71900 30.24300 26.67000
H  -2.23000 29.70200 25.13400
H  -0.77900 28.87000 25.77600
H  -3.74400 28.08300 30.96800
H  -0.90300 25.95500 33.37700
H  -1.91800 26.79500 34.52600
H  -0.58100 27.64700 33.62200


