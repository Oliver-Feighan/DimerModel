%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 53.03600 23.47300 43.74200
C  50.27300 25.68400 43.34900
C  50.89700 20.80500 43.54400
C  55.76900 21.39500 43.06400
C  55.19900 26.21400 43.34600
N  50.76700 23.35900 43.20300
C  49.85000 24.33200 43.33000
C  48.49000 23.71000 43.58900
C  48.66000 22.11500 43.58600
C  50.27200 22.11900 43.58300
C  48.04000 21.54600 44.88200
C  47.47500 24.17900 42.55900
C  47.98000 24.34500 41.10600
H  47.51400 23.28000 40.07600
N  53.32100 21.33300 43.43700
C  52.26200 20.45300 43.38400
C  52.84400 19.14800 43.34900
C  54.26800 19.24400 43.11100
C  54.55100 20.69600 43.23900
C  52.01300 17.88700 43.29400
C  55.22500 18.13300 42.67600
O  54.80000 17.01500 42.45700
C  56.63100 18.33100 42.48800
N  55.16000 23.75300 43.17000
C  56.04500 22.76900 43.10900
C  57.45600 23.24500 42.80400
C  57.16400 24.78400 42.53600
C  55.80000 24.97000 43.09600
C  58.49900 22.99600 44.00000
C  57.34100 25.12500 40.98500
C  58.47200 25.99800 40.52000
N  52.84200 25.47200 43.66600
C  53.82600 26.48100 43.63200
C  53.20700 27.81500 43.74300
C  51.85100 27.50000 43.65700
C  51.62900 26.08900 43.48500
C  53.87500 29.20800 44.01000
C  50.52800 28.13200 43.63200
O  50.14800 29.31100 43.72400
C  49.45600 27.00700 43.45200
C  48.88400 27.41400 42.13800
O  49.50700 27.53100 41.11400
O  47.60700 27.75600 42.28300
C  46.79800 28.12100 41.06800
H  50.22500 19.97800 43.78200
H  56.71400 20.86800 42.91800
H  55.94300 26.97900 43.11200
H  48.29500 24.16100 44.56200
H  48.23100 21.65000 42.69800
H  47.24600 20.85400 44.60200
H  47.57500 22.29200 45.52700
H  48.83900 20.99200 45.37400
H  47.06700 25.14300 42.86300
H  46.65900 23.46100 42.47200
H  49.05000 24.24200 40.92600
H  47.63400 25.33600 40.81200
H  50.94600 18.10600 43.24600
H  52.31600 17.25700 44.13100
H  52.35200 17.31500 42.43100
H  56.90400 19.01600 41.68500
H  57.10100 17.37500 42.25800
H  56.94000 18.76400 43.44000
H  57.86000 22.86500 41.86600
H  58.02800 25.24600 43.01400
H  58.00400 22.45500 44.80700
H  58.78500 23.97500 44.38500
H  59.43300 22.50300 43.72600
H  56.50400 25.71600 40.61300
H  57.40500 24.21500 40.39000
H  58.17800 27.02600 40.30800
H  58.72900 25.47900 39.59600
H  59.39500 25.89600 41.09100
H  54.34600 29.53100 43.08200
H  54.56000 29.07800 44.84800
H  53.08700 29.93100 44.22000
H  48.79600 27.08200 44.31700
H  46.90900 29.16800 40.78500
H  45.73600 27.93300 41.23100
H  47.15900 27.45300 40.28600
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


