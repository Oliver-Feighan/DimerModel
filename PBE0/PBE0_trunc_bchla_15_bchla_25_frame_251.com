%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg -2.51200 34.32000 26.71000
C  -3.70600 32.37000 29.49400
C  -1.09600 36.45200 28.85600
C  -2.23500 36.47200 24.15600
C  -4.45400 32.24200 24.65700
N  -2.47200 34.41100 28.88300
C  -2.84200 33.43600 29.82500
C  -2.27900 33.77800 31.25700
C  -1.51600 35.11400 31.00700
C  -1.64900 35.32300 29.47500
C  -1.96800 36.32500 31.78300
C  -1.39500 32.57400 31.88800
C  -1.81300 32.12300 33.29300
H  -0.66900 32.23700 34.36800
N  -1.63700 36.15800 26.50500
C  -1.11000 36.89500 27.48800
C  -0.61000 38.17800 26.99100
C  -1.04200 38.20400 25.62000
C  -1.69000 36.91400 25.36100
C  0.10900 39.14600 27.82900
C  -0.90500 39.45100 24.69100
O  -1.22300 39.41800 23.51400
C  -0.35500 40.75800 25.16600
N  -3.22300 34.33100 24.76300
C  -2.85500 35.30500 23.84100
C  -3.39300 34.95800 22.41400
C  -3.78700 33.44000 22.55100
C  -3.77600 33.28900 24.07200
C  -4.45800 35.91300 21.78300
C  -2.78800 32.46000 21.87700
C  -3.33300 31.25600 21.15100
N  -3.86100 32.67700 26.98500
C  -4.55900 31.95500 26.05900
C  -5.21600 30.76800 26.76200
C  -4.86600 30.86600 28.08700
C  -4.10200 32.07900 28.18800
C  -5.97600 29.69500 26.01600
C  -4.98700 30.30300 29.42000
O  -5.62100 29.26000 29.68400
C  -4.28200 31.35300 30.46900
C  -5.23900 32.03400 31.42500
O  -6.23700 32.64900 31.00000
O  -4.96700 31.73200 32.76000
C  -6.09800 32.02100 33.63400
H  -0.51900 37.13000 29.48800
H  -2.04400 37.20300 23.36800
H  -4.93400 31.57200 23.94000
H  -3.19000 33.89300 31.84400
H  -0.52800 34.94900 31.43900
H  -2.83400 36.13100 32.41600
H  -2.36000 36.92300 30.96100
H  -1.18900 36.74500 32.42000
H  -0.38600 32.98500 31.87000
H  -1.36700 31.74700 31.17900
H  -2.21200 31.11400 33.18900
H  -2.62300 32.76200 33.64700
H  -0.25500 40.16900 27.91900
H  1.01900 39.36100 27.27000
H  0.49900 38.75500 28.76900
H  0.66500 40.53600 25.48200
H  -0.91700 41.02300 26.06200
H  -0.46600 41.44500 24.32800
H  -2.50600 35.02500 21.78400
H  -4.81400 33.23900 22.24700
H  -4.50900 36.75000 22.47900
H  -5.42000 35.44200 21.58300
H  -4.08600 36.40600 20.88400
H  -2.24500 32.16900 22.77600
H  -2.06500 32.93000 21.21000
H  -2.75000 31.18100 20.23300
H  -4.37700 31.36900 20.85700
H  -3.03100 30.39100 21.74200
H  -5.97600 29.97600 24.96300
H  -7.01700 29.73400 26.33600
H  -5.54700 28.71900 26.24200
H  -3.48200 30.82800 30.99100
H  -6.59400 31.14900 34.05800
H  -6.86800 32.64200 33.17600
H  -5.69700 32.60500 34.46200


