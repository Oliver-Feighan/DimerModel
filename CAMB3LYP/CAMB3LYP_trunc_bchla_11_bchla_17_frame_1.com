%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 52.53200 24.27200 45.27800
C  49.69500 26.14500 43.97100
C  51.04500 21.43900 44.03200
C  55.55800 22.86200 45.08400
C  53.78600 27.15600 46.56500
N  50.65200 23.90000 43.94600
C  49.60000 24.82500 43.60600
C  48.41500 24.05300 43.01900
C  48.70000 22.56700 43.38700
C  50.22600 22.57000 43.84000
C  47.76900 22.00700 44.56400
C  48.04400 24.34300 41.59400
C  49.09500 25.06000 40.71400
H  48.96900 25.39700 39.21700
N  53.20200 22.33000 44.76700
C  52.42300 21.31600 44.25500
C  53.27600 20.09300 44.13100
C  54.57100 20.56600 44.32700
C  54.53600 21.98000 44.69300
C  52.68600 18.78600 43.73500
C  55.85700 19.70500 44.16100
O  55.72900 18.56500 43.75100
C  57.24900 20.21900 44.32400
N  54.40600 25.03400 45.62800
C  55.48900 24.15600 45.60300
C  56.74800 24.82300 46.15700
C  56.16900 26.07400 46.91200
C  54.72900 26.07700 46.39400
C  57.69100 23.98800 47.06100
C  56.92500 27.42300 46.68900
C  57.79700 27.94700 47.92000
N  51.87600 26.25900 45.36900
C  52.48600 27.32200 46.01900
C  51.58600 28.44900 46.02500
C  50.55800 28.02100 45.17400
C  50.79100 26.66600 44.83200
C  51.91200 29.79900 46.68700
C  49.30500 28.37700 44.58200
O  48.71800 29.46600 44.74300
C  48.58400 27.14400 43.87700
C  48.23600 27.64400 42.53200
O  48.98500 28.23500 41.74700
O  46.97500 27.21700 42.27700
C  46.32700 27.76200 41.04300
H  50.57900 20.54200 43.62000
H  56.53800 22.38300 45.02800
H  54.16100 28.03300 47.09600
H  47.65100 24.44400 43.69100
H  48.71700 21.90900 42.51800
H  47.15300 21.26400 44.05900
H  47.21300 22.83200 45.01000
H  48.44700 21.51600 45.26300
H  47.08700 24.86300 41.63300
H  47.84600 23.33500 41.22800
H  49.99900 24.46000 40.82300
H  49.29400 26.05200 41.12000
H  53.06600 18.57900 42.73500
H  51.60100 18.88400 43.75400
H  52.95100 18.03100 44.47600
H  57.60500 20.97300 43.62300
H  58.00300 19.45100 44.15400
H  57.23000 20.50700 45.37500
H  57.21200 25.04800 45.19600
H  56.15200 25.77400 47.96000
H  58.56700 23.60300 46.53800
H  57.16500 23.10200 47.41500
H  57.94800 24.60600 47.92200
H  56.28400 28.25100 46.38500
H  57.54500 27.23200 45.81300
H  58.71300 28.26700 47.42400
H  58.01300 27.19900 48.68300
H  57.20100 28.70700 48.42400
H  51.05000 30.36300 47.04400
H  52.35800 30.33200 45.84800
H  52.47600 29.64100 47.60600
H  47.71400 26.85200 44.46700
H  45.32600 27.34000 40.94700
H  46.88700 27.48700 40.14900
H  46.32300 28.83200 41.24600
Mg 29.46400 59.14500 41.03800
C  26.57800 57.38900 40.28900
C  31.34900 56.53100 39.78200
C  32.15000 61.18900 41.08600
C  27.35400 61.82000 41.94900
N  28.99200 57.29200 39.87100
C  27.77200 56.74800 39.89300
C  27.77900 55.38900 39.12700
C  29.32900 54.99700 39.21800
C  29.96700 56.34400 39.73200
C  29.58100 53.72700 40.06600
C  27.37400 55.59400 37.68500
C  26.46000 54.52700 37.10300
H  26.57800 54.02100 35.61600
N  31.53300 58.78800 40.65000
C  32.11000 57.71400 40.08400
C  33.51500 57.95800 40.03300
C  33.80600 59.24700 40.57400
C  32.47900 59.81600 40.84000
C  34.47700 56.89100 39.66100
C  35.16200 59.85800 40.59300
O  36.13000 59.21400 40.22400
C  35.39700 61.26400 41.13200
N  29.71300 61.28200 41.56500
C  30.94400 61.84600 41.41400
C  30.85900 63.33400 41.52000
C  29.29200 63.52200 41.79300
C  28.69900 62.13600 41.84900
C  31.79900 63.96200 42.56600
C  28.51800 64.48000 40.78300
C  27.64200 65.63900 41.31400
N  27.43700 59.54700 41.14000
C  26.74400 60.62700 41.58900
C  25.41200 60.28000 41.81900
C  25.25400 58.98400 41.34800
C  26.48900 58.61100 40.88000
C  24.32700 61.25500 42.36300
C  24.37200 57.92800 41.16800
O  23.15900 57.91800 41.34600
C  25.14900 56.75500 40.39700
C  25.25800 55.58200 41.18700
O  25.99300 55.31600 42.14100
O  24.33400 54.68200 40.66500
C  24.17500 53.35900 41.24400
H  31.95200 55.65500 39.53500
H  32.91900 61.96200 41.15300
H  26.80500 62.68200 42.33400
H  27.08300 54.71900 39.63100
H  29.81700 54.88800 38.25000
H  28.66900 53.37400 40.54700
H  30.34000 53.93800 40.82000
H  30.00600 52.97200 39.40500
H  28.16800 55.74800 36.95400
H  26.92700 56.58800 37.67000
H  25.44000 54.85100 37.30900
H  26.66400 53.69100 37.77200
H  34.81300 56.52500 40.63100
H  35.29700 57.28300 39.05900
H  34.02200 56.08300 39.08800
H  35.10600 61.93900 40.32700
H  36.47900 61.37300 41.20200
H  34.84400 61.42100 42.05800
H  31.14300 63.61000 40.50500
H  29.10900 64.04800 42.73000
H  32.59600 64.44300 41.99800
H  32.19100 63.23400 43.27700
H  31.26900 64.75100 43.09900
H  27.85100 63.86700 40.17700
H  29.18800 65.01800 40.11200
H  26.69000 65.11600 41.40400
H  27.54400 66.41700 40.55800
H  27.77900 66.00900 42.33000
H  24.06200 61.95200 41.56900
H  24.63400 61.90700 43.18100
H  23.42900 60.71200 42.65900
H  24.61600 56.61800 39.45600
H  23.13300 53.07400 41.10200
H  24.52900 53.25500 42.27000
H  24.72400 52.65900 40.61400


