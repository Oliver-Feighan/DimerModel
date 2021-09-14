%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 48.00100 15.36700 27.57000
C  45.75200 15.79400 30.49800
C  50.04400 17.47200 29.24900
C  49.70300 15.15900 24.90400
C  45.19600 13.91000 25.93200
N  47.85300 16.51500 29.59600
C  46.96000 16.37600 30.64500
C  47.41800 17.01300 31.96200
C  48.55800 17.91800 31.40100
C  48.86700 17.29700 29.98400
C  48.22700 19.37400 31.30400
C  47.82300 16.04400 33.12300
C  47.60700 16.50100 34.52700
H  48.68000 16.07500 35.52400
N  49.72200 16.17900 27.12000
C  50.48900 17.00200 27.97100
C  51.76400 17.20600 27.34600
C  51.69700 16.38300 26.13900
C  50.37000 15.89100 26.04500
C  52.97000 17.98900 27.98500
C  52.83500 16.13800 25.19600
O  52.69500 15.47100 24.14200
C  54.21900 16.61300 25.45000
N  47.52500 14.80200 25.63500
C  48.43700 14.77900 24.66100
C  47.90800 14.19700 23.31100
C  46.65000 13.42700 23.86100
C  46.39500 14.11500 25.23000
C  47.50200 15.31900 22.30200
C  46.93500 11.88100 24.07400
C  47.49500 11.04500 22.90400
N  45.91900 14.88800 28.08600
C  44.96000 14.24100 27.31500
C  43.72700 14.13900 28.06000
C  44.00400 14.69600 29.31500
C  45.32400 15.22300 29.27200
C  42.47900 13.44100 27.63500
C  43.42900 14.97200 30.62500
O  42.35300 14.72000 31.15100
C  44.65500 15.63100 31.44600
C  44.12500 16.85400 32.02400
O  43.69700 17.84400 31.41000
O  43.76300 16.58100 33.31200
C  43.14800 17.60200 34.09800
H  50.71300 18.03200 29.90600
H  50.25500 15.06200 23.96700
H  44.37700 13.44100 25.38200
H  46.60300 17.64800 32.31000
H  49.48300 17.75400 31.95300
H  48.08000 19.59600 30.24700
H  49.08100 19.96700 31.63100
H  47.34400 19.72300 31.83900
H  48.89700 15.88200 33.03600
H  47.33700 15.08800 32.93000
H  46.66100 16.03500 34.80200
H  47.32400 17.55200 34.58700
H  53.79100 17.37200 28.34900
H  52.63200 18.67300 28.76300
H  53.29200 18.60300 27.14300
H  54.25300 17.68500 25.25400
H  54.84800 16.07400 24.74100
H  54.46700 16.31100 26.46700
H  48.71600 13.57900 22.91800
H  45.74000 13.52600 23.27000
H  48.21500 15.26400 21.47900
H  47.43200 16.32600 22.71200
H  46.52600 14.96500 21.96700
H  45.99300 11.38800 24.31300
H  47.54700 11.81000 24.97300
H  48.43400 10.64000 23.28100
H  47.70000 11.78200 22.12700
H  46.83000 10.23200 22.61300
H  42.78200 12.85500 26.76700
H  41.77900 14.22600 27.34900
H  42.00600 12.85500 28.42300
H  44.98300 14.94300 32.22500
H  43.06700 18.54100 33.55000
H  43.82900 17.76200 34.93400
H  42.21500 17.25400 34.53900
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


