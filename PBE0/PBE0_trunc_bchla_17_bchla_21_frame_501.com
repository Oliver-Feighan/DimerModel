%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg 15.62100 51.82400 26.04800
C  17.47400 50.60700 28.62800
C  13.11000 52.75700 28.20700
C  13.81500 52.82700 23.36600
C  17.99400 50.42600 23.83300
N  15.23600 51.49000 28.16800
C  16.26900 51.27600 29.02600
C  15.94600 51.61800 30.44100
C  14.40600 52.07300 30.26100
C  14.21100 52.11900 28.78100
C  13.38100 51.07100 30.96200
C  16.88200 52.82900 30.90600
C  17.11300 53.00000 32.40900
H  16.27700 52.29800 33.50800
N  13.77400 52.76600 25.77300
C  12.87200 53.04400 26.81300
C  11.70500 53.64500 26.21000
C  11.90800 53.69900 24.81600
C  13.23300 53.12400 24.62900
C  10.53200 54.15900 27.07700
C  10.86200 54.19000 23.68400
O  11.18100 54.33400 22.49700
C  9.41300 54.55600 23.96000
N  15.93700 51.69300 23.86100
C  14.99700 52.16200 23.00600
C  15.27300 51.71300 21.57200
C  16.82600 51.21000 21.68000
C  17.02200 51.18500 23.20400
C  14.21700 50.69700 20.93700
C  17.94900 52.21900 21.11400
C  19.01900 51.70400 20.16900
N  17.42400 50.77200 26.12100
C  18.25400 50.26800 25.19100
C  19.45700 49.67400 25.85400
C  19.17000 49.86100 27.22700
C  17.86900 50.41600 27.30000
C  20.69200 49.07400 25.28100
C  19.74000 49.61300 28.49800
O  20.83000 49.08800 28.80500
C  18.66300 50.08200 29.47300
C  18.25000 48.96100 30.28600
O  17.83600 47.88700 29.97900
O  18.30900 49.39800 31.60600
C  17.88900 48.52200 32.76600
H  12.38500 53.27700 28.83600
H  13.18500 53.06500 22.50600
H  18.65800 49.80900 23.22500
H  16.01600 50.87800 31.23800
H  14.23400 53.06400 30.68100
H  12.50100 51.52900 31.41300
H  13.92900 50.42000 31.64300
H  12.99700 50.42600 30.17100
H  16.51300 53.78500 30.53400
H  17.84200 52.64500 30.42300
H  17.07200 54.07500 32.58900
H  18.14900 52.74600 32.63000
H  10.69700 54.07400 28.15100
H  9.65000 53.59900 26.76800
H  10.36400 55.20000 26.80000
H  9.47200 55.48000 24.53400
H  8.95500 53.91200 24.71100
H  8.82100 54.48200 23.04800
H  15.27500 52.60400 20.94300
H  16.96000 50.24900 21.18300
H  14.69800 49.77000 20.62400
H  13.70900 51.16400 20.09400
H  13.48400 50.38900 21.68300
H  18.44000 52.77000 21.91600
H  17.40700 53.00700 20.59200
H  19.97700 51.93800 20.63200
H  18.94100 52.12500 19.16700
H  18.97400 50.61500 20.12300
H  20.76700 49.12700 24.19500
H  20.76600 48.01700 25.53800
H  21.55200 49.59200 25.70500
H  19.01800 50.87900 30.12700
H  17.50200 47.57000 32.40200
H  17.07000 48.99900 33.30600
H  18.78800 48.39000 33.36800


