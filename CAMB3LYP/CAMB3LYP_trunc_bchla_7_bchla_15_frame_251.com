%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 25.94600 0.51000 29.41700
C  27.75500 0.04800 32.42200
C  23.08400 0.98300 31.28700
C  24.18800 0.65700 26.52400
C  28.83500 -0.18800 27.55800
N  25.56000 0.39300 31.59400
C  26.46800 0.41400 32.65000
C  25.73500 0.61600 34.03900
C  24.24300 0.62700 33.61200
C  24.25900 0.61600 32.06700
C  23.30400 -0.50900 34.16600
C  26.25400 1.80200 34.95100
C  26.76400 1.39400 36.30000
H  26.45300 2.38700 37.48900
N  23.95800 0.76900 29.03000
C  22.92100 1.01800 29.90600
C  21.67800 1.14300 29.16100
C  21.94400 1.05800 27.76500
C  23.43600 0.83800 27.70800
C  20.32000 1.49500 29.80100
C  21.00100 1.19500 26.57400
O  21.36800 1.18900 25.41700
C  19.44600 1.46900 26.69600
N  26.43900 0.16900 27.31600
C  25.54800 0.29300 26.34100
C  26.18900 0.06900 24.92500
C  27.73100 -0.01000 25.26000
C  27.69500 0.01500 26.81800
C  25.54500 -1.13400 24.17700
C  28.51900 1.10000 24.59300
C  29.85400 0.64900 24.00300
N  27.89800 -0.11200 29.79400
C  28.99000 -0.26300 28.96100
C  30.08800 -0.69100 29.77200
C  29.69800 -0.55100 31.12700
C  28.32100 -0.19300 31.11800
C  31.42300 -1.10700 29.28200
C  30.17600 -0.59300 32.52000
O  31.23000 -0.96300 33.02900
C  28.81000 -0.12800 33.44600
C  28.50200 -1.15100 34.39600
O  27.57900 -1.93000 34.29600
O  29.31100 -1.04800 35.50500
C  29.06500 -2.12600 36.49500
H  22.18200 1.26000 31.83700
H  23.66500 0.64600 25.56600
H  29.70300 -0.30600 26.90500
H  25.82900 -0.35200 34.53000
H  23.76600 1.56200 33.90600
H  22.35300 -0.03400 34.40900
H  23.77900 -0.95200 35.04100
H  23.10600 -1.31400 33.45800
H  25.42400 2.50500 35.02300
H  27.09400 2.25600 34.42500
H  27.85000 1.29700 36.26800
H  26.30600 0.41800 36.45800
H  20.08400 2.49100 29.42800
H  20.46000 1.63600 30.87300
H  19.49000 0.83300 29.55100
H  19.31200 2.38400 27.27400
H  18.87900 0.57700 26.96400
H  19.10200 1.53100 25.66400
H  26.08400 0.92600 24.26000
H  28.04800 -1.00900 24.96000
H  26.34500 -1.78700 23.82700
H  24.86300 -0.83400 23.38100
H  24.93000 -1.70700 24.86900
H  28.70300 1.93600 25.26700
H  27.98200 1.61200 23.79300
H  30.10400 -0.36600 24.31200
H  30.60900 1.32100 24.41200
H  29.90400 0.79600 22.92400
H  31.58400 -1.10700 28.20400
H  31.68500 -2.00300 29.84600
H  32.04200 -0.24000 29.51400
H  29.03700 0.80200 33.96800
H  29.24700 -3.09800 36.03600
H  28.05200 -1.96800 36.86500
H  29.74500 -2.08500 37.34500
Mg 46.81300 34.80600 27.83700
C  45.17600 32.79100 30.24000
C  46.71900 37.24900 30.02700
C  47.33800 36.80600 25.25700
C  46.11100 32.19400 25.53200
N  46.14800 34.85100 30.00200
C  45.62400 33.97200 30.77600
C  45.62400 34.45800 32.24000
C  45.67300 36.06300 32.08800
C  46.14500 36.09800 30.57500
C  44.31800 36.83800 32.33200
C  46.95000 33.94900 33.01300
C  46.90000 33.61900 34.54000
H  45.84300 34.41700 35.38900
N  47.21100 36.73800 27.68900
C  47.17400 37.59200 28.78100
C  47.54700 38.97300 28.31200
C  47.76600 38.81400 26.88000
C  47.45200 37.44200 26.47200
C  47.57300 40.14900 29.27000
C  48.24400 39.87200 25.91000
O  48.37500 39.70100 24.69000
C  48.62100 41.19700 26.43700
N  46.83600 34.49900 25.62400
C  47.04000 35.53300 24.82800
C  47.32500 35.03100 23.34300
C  46.87400 33.47700 23.48000
C  46.62500 33.36900 24.98200
C  46.58500 35.85900 22.24100
C  47.86700 32.37700 22.98000
C  47.46400 31.59700 21.75500
N  46.03500 32.76500 27.86100
C  45.83400 31.87700 26.82300
C  45.15900 30.73800 27.28800
C  44.84600 31.07700 28.61800
C  45.34900 32.31400 28.94500
C  44.79200 29.55000 26.55100
C  44.25600 30.51400 29.77000
O  43.77000 29.39300 30.08200
C  44.20600 31.72700 30.80400
C  44.55000 31.25500 32.16900
O  45.61400 30.77700 32.55700
O  43.43700 31.56400 33.01300
C  43.54100 31.17800 34.42600
H  46.79200 37.98600 30.83000
H  47.57400 37.46800 24.42200
H  45.97100 31.36600 24.83300
H  44.62300 34.28200 32.63300
H  46.32300 36.45800 32.86800
H  43.46200 36.16400 32.35700
H  44.09900 37.53000 31.51900
H  44.43400 37.33900 33.29300
H  47.78200 34.63400 32.84900
H  47.31500 33.05900 32.50000
H  47.84900 33.77500 35.05400
H  46.48500 32.61400 34.62600
H  48.46100 40.76800 29.13700
H  47.51700 39.85100 30.31700
H  46.63300 40.69500 29.19400
H  47.78800 41.66200 26.96500
H  48.94200 41.84200 25.61900
H  49.47200 41.20900 27.11800
H  48.38800 35.01100 23.10300
H  45.86800 33.30600 23.09400
H  45.91300 36.58500 22.69900
H  46.15500 35.05600 21.64200
H  47.35600 36.29200 21.60500
H  47.98100 31.65400 23.78800
H  48.84800 32.82300 22.81500
H  46.58600 32.13300 21.39700
H  47.14500 30.60700 22.08300
H  48.21600 31.60800 20.96700
H  45.69700 29.19600 26.05800
H  44.17700 29.70800 25.66500
H  44.42700 28.78900 27.24100
H  43.19600 32.13800 30.79400
H  42.60700 31.38000 34.95100
H  44.23500 31.85900 34.91900
H  43.73000 30.10700 34.49100


