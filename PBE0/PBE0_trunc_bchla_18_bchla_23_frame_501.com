%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg -10.08900 40.79200 42.94700
C  -8.61600 37.87300 41.73900
C  -7.78200 42.64100 41.17400
C  -11.69100 43.40600 44.11700
C  -13.01100 38.80800 43.84600
N  -8.29800 40.23100 41.61500
C  -7.89900 38.98100 41.23300
C  -6.67800 39.04800 40.38600
C  -6.44000 40.59600 40.17200
C  -7.57700 41.25700 41.04500
C  -4.98600 41.06900 40.42300
C  -6.83200 38.29900 39.04000
C  -5.67600 37.32300 38.82800
H  -4.88500 37.42500 37.48900
N  -9.75400 42.76600 42.70300
C  -8.72900 43.36600 41.97300
C  -8.87200 44.79700 42.20800
C  -10.00600 45.05500 43.11000
C  -10.52500 43.68700 43.42700
C  -8.10300 45.78800 41.48800
C  -10.49500 46.41500 43.39200
O  -9.97600 47.42800 42.97200
C  -11.74500 46.70800 44.14600
N  -12.01700 41.03600 43.77600
C  -12.30900 42.20100 44.43900
C  -13.68100 42.11000 45.07300
C  -14.24200 40.84200 44.54700
C  -13.08300 40.19800 43.96600
C  -13.70800 42.29100 46.64600
C  -15.26700 41.05700 43.38100
C  -16.32200 40.05900 43.24400
N  -10.72900 38.74300 42.86600
C  -11.88000 38.08300 43.31600
C  -11.64700 36.67300 43.29400
C  -10.42100 36.52600 42.65600
C  -9.90000 37.83800 42.36900
C  -12.60800 35.66600 43.80900
C  -9.50700 35.58200 42.14800
O  -9.55000 34.37200 42.15700
C  -8.30600 36.38400 41.59900
C  -7.07900 35.93500 42.27600
O  -6.65100 36.36600 43.36000
O  -6.55600 34.85100 41.60900
C  -5.28400 34.28000 42.02100
H  -6.98400 43.15900 40.63900
H  -12.35200 44.18100 44.51300
H  -13.89800 38.25100 44.15500
H  -5.86700 38.79100 41.06800
H  -6.60400 40.92400 39.14600
H  -4.42500 40.13600 40.35700
H  -4.81000 41.64900 41.32900
H  -4.66100 41.71400 39.60700
H  -6.76000 38.94200 38.16300
H  -7.79900 37.79800 39.08100
H  -6.24400 36.39600 38.75100
H  -5.06200 37.35000 39.72700
H  -7.39900 45.26300 40.84200
H  -7.46100 46.42800 42.09400
H  -8.77100 46.37100 40.85300
H  -11.59200 46.04300 44.99600
H  -12.59700 46.39100 43.54400
H  -11.78600 47.75300 44.45400
H  -14.19200 42.98500 44.67300
H  -14.60800 40.10900 45.26600
H  -12.74000 42.63900 47.00700
H  -14.11000 41.44200 47.20000
H  -14.30700 43.14200 46.96900
H  -14.72400 41.14600 42.44000
H  -15.78500 42.00900 43.49600
H  -17.15500 40.34300 43.88800
H  -15.96800 39.05000 43.45800
H  -16.80100 40.16400 42.27000
H  -13.09400 36.01500 44.72000
H  -12.17300 34.67900 43.96700
H  -13.39900 35.66300 43.05900
H  -8.06700 36.13500 40.56500
H  -5.23500 33.70700 42.94700
H  -4.52300 35.06000 42.03400
H  -5.07200 33.53500 41.25500


