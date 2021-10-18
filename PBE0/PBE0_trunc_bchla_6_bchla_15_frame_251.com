%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 17.10100 -2.07700 27.77000
C  16.45500 -0.19100 30.58100
C  19.23600 -3.94700 29.68100
C  17.96100 -3.74600 24.96300
C  15.37700 0.23000 25.82200
N  17.83400 -1.91300 29.91400
C  17.31200 -1.24300 30.86800
C  17.86400 -1.55200 32.35400
C  18.99500 -2.58100 31.95700
C  18.76000 -2.80300 30.38800
C  20.50900 -2.04800 32.15900
C  16.79200 -2.13400 33.37400
C  16.96200 -1.78000 34.89400
H  18.25700 -1.14000 35.29400
N  18.31900 -3.63600 27.38700
C  19.10500 -4.38300 28.35000
C  19.74900 -5.47800 27.63200
C  19.48000 -5.36700 26.24500
C  18.51900 -4.24500 26.14900
C  20.50800 -6.46400 28.44500
C  19.91900 -6.24400 25.04900
O  19.71600 -6.00400 23.89600
C  20.84100 -7.42500 25.41400
N  16.84100 -1.78700 25.63700
C  17.26600 -2.56700 24.69400
C  17.04700 -1.94200 23.28700
C  16.12300 -0.73100 23.56000
C  16.07300 -0.74300 25.10600
C  18.39300 -1.56900 22.63000
C  14.74200 -0.86000 22.84700
C  14.78100 -0.68700 21.32000
N  16.09100 -0.30700 28.09200
C  15.35600 0.45700 27.20300
C  14.73300 1.49900 27.97600
C  15.14100 1.35600 29.33600
C  15.94900 0.18700 29.29900
C  13.91400 2.52600 27.44500
C  15.02300 1.80700 30.68600
O  14.43000 2.83800 31.12300
C  15.98100 0.86900 31.51800
C  16.98700 1.71200 32.12800
O  17.99700 2.11700 31.55300
O  16.51500 2.07400 33.36400
C  17.31500 3.14200 34.00000
H  19.94300 -4.54500 30.25900
H  18.26400 -4.34800 24.10300
H  14.78900 0.88500 25.17500
H  18.31000 -0.60500 32.65900
H  18.90200 -3.48800 32.55400
H  20.47800 -0.95900 32.20400
H  21.12500 -2.36600 31.31800
H  21.03200 -2.40200 33.04800
H  16.90400 -3.21300 33.27100
H  15.74800 -2.07700 33.06400
H  16.83100 -2.63200 35.56100
H  16.22600 -1.02300 35.16500
H  21.52400 -6.11200 28.62600
H  20.42400 -7.48300 28.06700
H  20.00400 -6.52900 29.40900
H  20.33800 -8.26100 25.90100
H  21.56700 -6.93500 26.06200
H  21.42700 -7.84100 24.59500
H  16.51200 -2.73300 22.76200
H  16.63000 0.18800 23.26600
H  19.22100 -1.38900 23.31600
H  18.44300 -0.64500 22.05300
H  18.63700 -2.41000 21.98100
H  14.08400 -0.07400 23.21500
H  14.28900 -1.82100 23.09300
H  15.82700 -0.60000 21.02600
H  14.42000 0.31100 21.07100
H  14.22200 -1.38300 20.69500
H  13.25200 2.01200 26.74800
H  14.42200 3.30300 26.87400
H  13.44500 3.05700 28.27400
H  15.32300 0.48100 32.29600
H  17.17000 4.03200 33.38700
H  18.34600 2.83300 34.16900
H  16.87700 3.33200 34.98000
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


