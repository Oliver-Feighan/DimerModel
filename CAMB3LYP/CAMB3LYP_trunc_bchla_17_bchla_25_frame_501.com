%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 29.51000 58.90800 41.21400
C  26.38900 57.48600 40.17300
C  31.07800 56.05600 39.82600
C  32.44700 60.51200 41.46200
C  27.84100 61.85400 41.88400
N  28.82700 56.95900 40.12100
C  27.51800 56.70000 39.82500
C  27.43600 55.53200 38.85500
C  28.91100 54.97600 39.03000
C  29.71500 56.02800 39.69400
C  29.08100 53.62400 39.80100
C  26.96300 56.03100 37.45200
C  25.85600 55.22200 36.65900
H  26.19700 54.60200 35.29200
N  31.43900 58.34900 40.91600
C  31.91100 57.13400 40.35300
C  33.36200 57.11100 40.33800
C  33.74700 58.43000 40.81400
C  32.51100 59.13800 41.09000
C  34.11400 55.98100 39.76900
C  35.21500 58.86100 40.97700
O  36.08800 58.03800 40.57400
C  35.69500 60.19600 41.58500
N  30.14600 61.00100 41.36300
C  31.38600 61.37000 41.62100
C  31.40400 62.84700 41.84900
C  29.90900 63.22000 42.03600
C  29.21400 61.95000 41.65200
C  32.31600 63.46200 42.99500
C  29.40100 64.48800 41.17400
C  28.93400 65.75600 41.93000
N  27.56300 59.51200 41.26600
C  27.03200 60.71300 41.55100
C  25.61900 60.69600 41.44600
C  25.33000 59.42300 40.94900
C  26.56000 58.73200 40.80800
C  24.61100 61.72700 41.76100
C  24.21900 58.63400 40.51400
O  23.05400 58.95300 40.43300
C  24.92200 57.26000 39.91000
C  24.42400 56.08500 40.63500
O  24.27500 56.00600 41.86100
O  24.19400 55.10700 39.69100
C  23.68300 53.83100 40.18900
H  31.62000 55.20100 39.41900
H  33.35900 61.11200 41.46000
H  27.25300 62.67600 42.29800
H  26.78700 54.79200 39.32400
H  29.19500 54.98000 37.97700
H  29.75800 53.81000 40.63600
H  29.53300 52.82700 39.21100
H  28.17800 53.23200 40.27000
H  27.91200 55.97700 36.91900
H  26.65200 57.07400 37.39700
H  25.02500 55.92000 36.55900
H  25.42300 54.48900 37.33900
H  33.71100 55.81900 38.76900
H  33.89700 55.07100 40.32800
H  35.20300 56.02600 39.79800
H  35.47200 60.98900 40.87100
H  36.68900 59.96900 41.96900
H  35.05300 60.30300 42.46000
H  31.80800 63.29500 40.94200
H  29.82800 63.26800 43.12200
H  32.55800 64.46300 42.63700
H  33.18900 62.82200 43.12200
H  31.69200 63.52000 43.88700
H  28.59500 64.07600 40.56700
H  30.17400 64.71200 40.43900
H  27.96900 66.12300 41.58300
H  29.69100 66.53100 41.80700
H  28.89100 65.63000 43.01200
H  23.61700 61.29700 41.88000
H  24.71500 62.31200 40.84700
H  25.04500 62.40800 42.49300
H  24.83500 57.21500 38.82500
H  22.80800 53.96700 40.82500
H  24.47100 53.30400 40.72900
H  23.37200 53.22500 39.33800
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


