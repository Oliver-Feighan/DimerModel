%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 40.80400 41.48400 27.15800
C  39.93800 43.96600 29.51800
C  41.67100 39.43500 29.59200
C  42.41200 39.59100 24.81200
C  40.80000 44.20800 24.65700
N  40.84000 41.74700 29.31900
C  40.23400 42.70800 30.04600
C  40.11400 42.23500 31.50800
C  40.62500 40.83000 31.52200
C  41.12800 40.63200 30.06000
C  41.62900 40.55000 32.71400
C  38.74000 42.40700 32.23800
C  38.20000 41.35400 33.29800
H  37.58800 41.87900 34.60500
N  41.67300 39.62800 27.19200
C  41.94100 38.92200 28.34300
C  42.50600 37.66000 27.95800
C  42.79100 37.65600 26.53700
C  42.34000 39.00600 26.13300
C  42.84100 36.69600 29.02300
C  43.34600 36.53900 25.61900
O  43.54900 36.75600 24.47000
C  43.39500 35.12400 26.09600
N  41.54800 41.90500 25.04000
C  42.04300 40.86300 24.33500
C  42.34300 41.24800 22.85100
C  41.79900 42.70700 22.76800
C  41.36100 43.00200 24.22400
C  43.78600 41.15300 22.43200
C  40.63300 43.00600 21.76200
C  39.24100 42.61000 22.12200
N  40.45200 43.63600 26.97900
C  40.44000 44.51600 25.93500
C  39.94600 45.79200 26.39200
C  39.72700 45.66500 27.78900
C  40.09700 44.34100 28.10100
C  39.86300 47.06500 25.62600
C  39.33400 46.33800 28.99400
O  38.92600 47.46900 29.23200
C  39.40400 45.25700 30.16400
C  40.16300 45.75100 31.38400
O  41.33900 46.00900 31.38400
O  39.30000 45.80400 32.46700
C  39.93400 46.30500 33.68700
H  41.79100 38.62500 30.31500
H  42.89000 38.90600 24.10900
H  40.69600 44.98800 23.90000
H  40.78700 42.96300 31.96100
H  39.78000 40.17800 31.74500
H  41.64600 41.41100 33.38200
H  42.57700 40.23500 32.27800
H  41.25100 39.76900 33.37400
H  37.98000 42.45900 31.45800
H  38.66600 43.37100 32.74200
H  39.06500 40.72800 33.51800
H  37.46000 40.71500 32.81600
H  43.75800 36.20700 28.69500
H  42.05100 35.94600 29.03300
H  42.98100 37.00500 30.05900
H  44.39900 35.05700 26.51300
H  43.13900 34.46100 25.26900
H  42.74200 34.89700 26.93900
H  41.71500 40.69400 22.15300
H  42.49200 43.53200 22.60600
H  44.38600 40.63600 23.18200
H  44.24800 42.13400 22.31800
H  43.78100 40.59500 21.49600
H  40.84900 42.53700 20.80200
H  40.70100 44.05800 21.48400
H  39.26600 41.80800 22.86000
H  38.56300 42.19300 21.37700
H  38.76600 43.53300 22.45600
H  40.81300 47.42600 25.23300
H  39.70300 47.84000 26.37600
H  39.09100 47.02600 24.85700
H  38.36800 45.04000 30.42500
H  39.84500 45.31900 34.14200
H  39.21300 47.03100 34.06200
H  40.94500 46.67900 33.52400
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


