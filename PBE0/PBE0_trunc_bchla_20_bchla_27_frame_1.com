%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 7.14200 57.54600 41.51700
C  6.28900 54.20700 41.28700
C  10.28800 56.72100 40.38900
C  7.81200 60.89700 41.16900
C  3.94100 58.30400 42.53200
N  8.10000 55.73200 40.72500
C  7.55000 54.46800 40.80800
C  8.51900 53.39300 40.32300
C  9.88500 54.17400 40.37000
C  9.39700 55.61300 40.45200
C  10.94100 53.73800 41.48400
C  8.18100 52.71900 38.95500
C  8.31100 51.21100 38.91500
H  9.48100 50.75900 38.05500
N  8.80400 58.60500 40.85800
C  9.98200 58.10200 40.47200
C  10.93200 59.16600 40.19700
C  10.33500 60.40200 40.38100
C  8.90800 60.02900 40.77200
C  12.41400 58.88100 39.65800
C  11.01900 61.80700 40.27800
O  12.22600 61.86700 39.90900
C  10.28000 63.06100 40.55000
N  5.92900 59.42300 41.58100
C  6.46000 60.62400 41.57700
C  5.53200 61.67900 42.14200
C  4.20400 60.94400 42.52900
C  4.70000 59.43900 42.28400
C  6.08200 62.75400 43.13300
C  2.88100 61.24800 41.76800
C  1.64800 60.78800 42.51800
N  5.41300 56.46200 41.84200
C  4.16700 56.89800 42.34900
C  3.34900 55.73600 42.56600
C  4.13900 54.62000 42.19700
C  5.33900 55.15600 41.73800
C  2.01600 55.69000 43.29800
C  4.11700 53.24300 42.14700
O  3.28800 52.39200 42.32200
C  5.50100 52.88200 41.39400
C  5.07100 52.59100 39.95700
O  4.83200 53.42500 39.10500
O  4.99200 51.21300 39.83600
C  4.48800 50.71300 38.56100
H  11.34700 56.50800 40.23000
H  8.09700 61.94400 41.04900
H  2.94700 58.37800 42.97800
H  8.57300 52.63000 41.09900
H  10.29300 54.14500 39.35900
H  11.53300 54.64700 41.59000
H  11.67500 53.00700 41.14400
H  10.53800 53.39000 42.43500
H  8.89800 53.17500 38.27200
H  7.26400 53.03800 38.46000
H  7.37100 50.81600 38.52900
H  8.48000 50.78000 39.90200
H  12.65000 57.83100 39.48000
H  13.15600 59.27900 40.35100
H  12.59000 59.28000 38.66000
H  9.90800 63.04900 41.57400
H  9.55300 63.27700 39.76700
H  11.06400 63.81300 40.63000
H  5.09100 62.35700 41.41100
H  4.06100 60.88100 43.60800
H  7.16900 62.67900 43.08000
H  5.65900 62.68700 44.13600
H  5.85900 63.78200 42.84700
H  2.99200 60.79700 40.78200
H  2.91000 62.32500 41.60300
H  0.83700 61.47100 42.26400
H  1.73700 60.78400 43.60500
H  1.32700 59.82700 42.11600
H  1.67100 54.68100 43.52400
H  1.28900 56.22800 42.69100
H  2.06800 56.22200 44.24800
H  6.04900 52.10100 41.92200
H  3.44600 51.01000 38.44100
H  4.61600 49.63500 38.46900
H  5.02300 51.19800 37.74500
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


