%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 16.41100 52.36700 25.56600
C  18.01800 50.88800 28.30000
C  13.92400 53.47100 27.63400
C  14.78500 53.43800 22.89200
C  18.97700 50.98900 23.49400
N  16.01700 52.13900 27.73500
C  16.81400 51.56600 28.70100
C  16.31600 51.93200 30.08700
C  14.87100 52.50500 29.80500
C  15.00600 52.79900 28.27200
C  13.75200 51.51700 30.10900
C  17.38000 52.82700 30.82100
C  16.86600 53.42500 32.16000
H  15.83800 52.61100 32.91800
N  14.67300 53.30500 25.31300
C  13.70900 53.59600 26.27000
C  12.51400 54.21500 25.67000
C  12.80900 54.23100 24.27200
C  14.13200 53.63200 24.11700
C  11.36800 54.71200 26.46900
C  11.91000 54.75300 23.18100
O  12.26900 54.77100 21.97100
C  10.56700 55.35400 23.49600
N  16.86500 52.27500 23.47100
C  15.99700 52.81000 22.61200
C  16.44300 52.50800 21.11400
C  17.81800 51.65600 21.31000
C  17.89300 51.61900 22.83700
C  15.27500 51.72000 20.39100
C  19.11700 52.36100 20.72800
C  19.40800 52.10700 19.22100
N  18.08500 51.15100 25.76500
C  19.11600 50.75400 24.87300
C  20.19800 50.00100 25.64600
C  19.78200 49.96200 26.97300
C  18.52400 50.75100 27.00100
C  21.48500 49.44800 25.11900
C  20.10700 49.55500 28.32000
O  21.06000 48.92400 28.78800
C  19.05800 50.30100 29.28500
C  18.45300 49.37700 30.31700
O  17.48000 48.69900 30.12500
O  19.02600 49.59100 31.59100
C  18.26600 49.12100 32.75900
H  13.19300 53.70900 28.41000
H  14.23600 53.78000 22.01200
H  19.73200 50.45400 22.91400
H  16.24800 51.01000 30.66300
H  14.63800 53.45400 30.28900
H  13.10300 52.02100 30.82500
H  14.20000 50.58000 30.43900
H  13.06200 51.39100 29.27500
H  17.73700 53.44800 30.00000
H  18.21600 52.16100 31.03200
H  16.46200 54.39100 31.85700
H  17.72400 53.56200 32.81800
H  10.55700 54.06200 26.14100
H  11.17100 55.77400 26.32500
H  11.38200 54.60100 27.55400
H  10.71400 56.25300 24.09300
H  10.12800 54.49800 24.01000
H  9.92000 55.56900 22.64600
H  16.62800 53.52500 20.76600
H  17.77900 50.58400 21.11200
H  14.41400 51.50600 21.02400
H  15.64500 50.76000 20.03000
H  14.90700 52.39400 19.61800
H  19.98300 52.23400 21.37700
H  18.92400 53.43100 20.79800
H  19.48300 53.08600 18.74900
H  18.62700 51.45900 18.82300
H  20.36100 51.60000 19.06800
H  22.27100 50.10500 25.49200
H  21.46400 49.41600 24.02900
H  21.78400 48.45300 25.44800
H  19.61000 51.11500 29.75400
H  18.44700 48.04700 32.79300
H  17.22900 49.45300 32.70500
H  18.68800 49.63700 33.62100
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


