%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
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
Mg 25.49300 50.21300 27.08700
C  23.37900 51.54000 29.65700
C  27.62800 49.17200 29.69500
C  27.85200 49.59100 24.90700
C  23.46200 51.62200 24.75200
N  25.44300 50.31400 29.39800
C  24.41900 50.80900 30.13400
C  24.66700 50.47600 31.59600
C  26.00400 49.68800 31.68900
C  26.38100 49.68400 30.19500
C  27.10100 50.37400 32.68300
C  23.49400 49.61900 32.23900
C  22.96400 50.09200 33.59300
H  23.92100 50.77900 34.64800
N  27.50700 49.52400 27.28400
C  28.12900 49.15000 28.38200
C  29.45200 48.73800 28.00600
C  29.55900 48.68700 26.59300
C  28.31600 49.30500 26.21200
C  30.44200 48.13300 29.00000
C  30.57100 48.09600 25.62400
O  30.49100 48.23300 24.41700
C  31.66500 47.26500 26.23600
N  25.76300 50.88600 25.17100
C  26.76300 50.29700 24.39500
C  26.51800 50.44800 22.86900
C  25.02700 50.95500 22.90400
C  24.70800 51.17300 24.38300
C  27.54000 51.40800 22.24300
C  23.97400 49.88700 22.34200
C  23.37400 50.12000 20.99600
N  23.75000 51.22900 27.14200
C  23.04900 51.77400 26.08800
C  21.90700 52.57500 26.52100
C  21.98600 52.51200 27.95800
C  23.17000 51.80300 28.27000
C  20.85200 53.08800 25.60900
C  21.31100 52.84500 29.20300
O  20.22900 53.44600 29.33300
C  22.22800 52.23700 30.41000
C  22.63900 53.39700 31.18500
O  23.23000 54.34500 30.71600
O  22.29300 53.27500 32.47200
C  22.54800 54.35700 33.47400
H  28.29000 48.76300 30.46100
H  28.36300 49.05600 24.10300
H  22.79700 51.98300 23.96500
H  24.78800 51.42300 32.12200
H  25.86900 48.70800 32.14700
H  26.61300 51.02100 33.41200
H  27.77100 50.92200 32.02000
H  27.67200 49.59700 33.19100
H  23.69200 48.55100 32.31900
H  22.62100 49.66400 31.58800
H  22.46400 49.28700 34.13200
H  22.16000 50.77100 33.30900
H  30.45600 47.04300 29.00900
H  30.32700 48.56400 29.99500
H  31.42700 48.48300 28.69200
H  31.35000 46.60400 27.04400
H  32.50000 47.83500 26.64100
H  31.97600 46.50800 25.51600
H  26.62800 49.45900 22.42500
H  24.97100 51.86600 22.30800
H  28.36600 51.68000 22.90000
H  27.08200 52.35600 21.96300
H  27.94500 50.95800 21.33600
H  23.09600 49.73600 22.97000
H  24.37200 48.87300 22.33000
H  23.66900 49.32700 20.31000
H  23.80100 51.07900 20.70200
H  22.28600 50.13300 21.06000
H  19.87400 53.08000 26.08900
H  20.94000 52.51200 24.68800
H  21.11300 54.12300 25.38400
H  21.65900 51.56400 31.05000
H  22.17700 55.30500 33.08400
H  23.60800 54.48100 33.69700
H  22.16500 54.13200 34.46900


