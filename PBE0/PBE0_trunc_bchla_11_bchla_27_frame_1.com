%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg -5.21000 24.54100 27.19500
C  -3.54100 26.38400 29.69300
C  -6.38000 22.49900 29.72600
C  -6.78300 22.59100 24.87600
C  -3.73200 26.39300 24.81700
N  -5.00900 24.47500 29.48500
C  -4.39300 25.40300 30.31600
C  -4.77800 25.15700 31.76100
C  -5.38000 23.70500 31.74000
C  -5.61200 23.56800 30.25200
C  -4.48700 22.58300 32.40100
C  -5.82800 26.14800 32.28500
C  -5.66000 26.57200 33.74800
H  -5.92800 25.52000 34.86000
N  -6.48000 22.85300 27.28100
C  -6.76500 22.19100 28.42600
C  -7.52200 21.02600 28.10100
C  -7.81700 21.04500 26.69400
C  -6.98500 22.14200 26.23200
C  -7.94700 20.03100 29.14500
C  -8.63000 20.10000 25.80900
O  -8.90700 20.38500 24.62300
C  -9.08700 18.72400 26.28000
N  -5.32000 24.55800 25.17000
C  -5.97300 23.61500 24.38200
C  -5.69900 23.91000 22.94300
C  -4.77200 25.17400 22.83700
C  -4.60800 25.45200 24.34000
C  -5.24900 22.67900 22.10000
C  -5.37300 26.38000 22.08300
C  -4.94400 26.45700 20.60500
N  -3.89500 26.01800 27.18100
C  -3.38100 26.77900 26.13200
C  -2.50100 27.90900 26.61000
C  -2.55000 27.71500 28.00000
C  -3.42600 26.62300 28.30100
C  -1.78700 28.93500 25.74900
C  -2.03300 28.27400 29.24000
O  -1.24900 29.20200 29.44300
C  -2.78100 27.47500 30.39000
C  -1.66500 26.98400 31.33300
O  -0.58100 26.54900 31.00300
O  -2.03300 27.27600 32.64500
C  -1.03300 26.85400 33.62600
H  -6.79600 21.92300 30.55500
H  -7.24600 22.04300 24.05300
H  -3.20100 26.86300 23.98700
H  -3.94200 25.15700 32.46000
H  -6.32400 23.59600 32.27400
H  -3.58700 22.93500 32.90600
H  -4.09000 21.81900 31.73300
H  -5.03500 22.07100 33.19200
H  -6.83700 25.75100 32.17600
H  -5.77900 27.08800 31.73500
H  -6.23400 27.46900 33.98000
H  -4.65000 26.93300 33.94100
H  -7.01100 19.85400 29.67500
H  -8.29300 19.05500 28.80500
H  -8.67800 20.40400 29.86200
H  -8.25700 18.16700 26.71500
H  -9.45700 18.28200 25.35500
H  -9.85300 18.80700 27.05100
H  -6.66900 24.06000 22.46800
H  -3.77200 24.98400 22.44500
H  -5.16100 21.82900 22.77700
H  -4.26700 22.88500 21.67500
H  -6.05700 22.54500 21.38000
H  -5.16500 27.26900 22.67800
H  -6.45400 26.27700 21.98600
H  -4.65800 25.46600 20.25200
H  -4.05700 27.05100 20.38400
H  -5.81000 26.77000 20.02100
H  -2.54700 29.61400 25.36200
H  -1.26100 28.53800 24.88100
H  -1.10500 29.45700 26.42000
H  -3.50500 28.04800 30.96900
H  -1.16000 25.83800 34.00100
H  -1.12400 27.59400 34.42100
H  0.01800 26.83100 33.33800


