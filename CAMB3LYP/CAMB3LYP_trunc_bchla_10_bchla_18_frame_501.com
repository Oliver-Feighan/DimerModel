%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 40.46600 8.46800 29.02600
C  42.13100 9.82400 31.82300
C  38.27200 6.93100 31.17900
C  38.87200 7.07400 26.41900
C  43.04200 9.56800 26.92900
N  40.22800 8.55200 31.19800
C  40.93100 9.11600 32.17600
C  40.40800 8.90400 33.61200
C  39.36400 7.73500 33.26600
C  39.31400 7.70200 31.75800
C  39.74100 6.39300 33.84500
C  39.69500 10.21300 34.26600
C  39.76000 10.23100 35.79600
H  40.52600 11.44400 36.38400
N  38.67500 7.29000 28.78400
C  37.94400 6.71700 29.85400
C  36.72400 6.18800 29.28700
C  36.85500 6.23500 27.90200
C  38.16400 6.90300 27.61900
C  35.64600 5.49400 30.09800
C  35.92100 5.81300 26.79100
O  36.21800 6.02700 25.65200
C  34.58600 5.19600 27.06600
N  40.97200 8.15700 27.00800
C  40.11100 7.61300 26.10000
C  40.78500 7.77300 24.70900
C  41.70800 8.97300 24.82300
C  41.98000 8.92500 26.34500
C  41.61300 6.52600 24.40200
C  41.12800 10.25100 24.34100
C  39.73300 10.63200 24.77100
N  42.22500 9.51500 29.25500
C  43.18900 9.89200 28.33900
C  44.24200 10.64200 29.04100
C  43.84000 10.65000 30.38800
C  42.65100 9.97700 30.46100
C  45.46200 11.19800 28.44400
C  44.09900 11.19900 31.69200
O  45.01500 11.97300 32.09800
C  42.93800 10.67900 32.66500
C  43.61600 9.87100 33.74100
O  44.43500 9.03300 33.50500
O  43.17500 10.12500 35.06800
C  43.60100 9.11400 36.08300
H  37.59700 6.39300 31.84700
H  38.44200 6.55400 25.56000
H  43.78300 10.08700 26.31700
H  41.23500 8.65200 34.27700
H  38.42900 8.08400 33.70600
H  39.16600 6.28900 34.76500
H  40.81800 6.26800 33.95900
H  39.27800 5.67800 33.16500
H  38.64200 10.24000 33.98400
H  40.12300 11.12800 33.85600
H  40.22100 9.26500 36.00300
H  38.77200 10.15700 36.25000
H  35.82300 5.83000 31.11900
H  35.93800 4.44600 30.03700
H  34.64100 5.77100 29.77900
H  34.06100 6.07300 27.44500
H  34.59100 4.42700 27.83800
H  34.14900 4.81100 26.14400
H  40.09900 7.86800 23.86800
H  42.63400 8.90500 24.25300
H  42.69600 6.55200 24.52300
H  41.40100 6.24300 23.37100
H  41.23400 5.71200 25.02100
H  41.17500 10.33500 23.25500
H  41.80700 11.00800 24.73300
H  39.85400 11.49300 25.42900
H  39.05700 9.95800 25.29800
H  39.19800 10.99600 23.89400
H  46.11400 10.36800 28.17200
H  45.85400 11.85300 29.22200
H  45.24700 11.72600 27.51400
H  42.43100 11.58700 32.99200
H  44.64200 9.27500 36.36300
H  43.48400 8.08500 35.74000
H  43.10400 9.17100 37.05100
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


