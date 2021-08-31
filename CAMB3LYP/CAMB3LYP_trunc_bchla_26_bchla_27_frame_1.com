%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg -9.12600 18.73000 43.07000
C  -5.62300 18.44500 42.86700
C  -8.71100 22.08300 42.24600
C  -12.43100 19.03300 42.79300
C  -9.30300 15.26600 42.95700
N  -7.32000 20.02500 42.45300
C  -5.96700 19.76600 42.55400
C  -5.07000 20.98600 42.54900
C  -6.06600 22.15700 42.37300
C  -7.44800 21.42600 42.33500
C  -5.98600 23.25800 43.50700
C  -3.93600 21.03100 41.43000
C  -4.22600 20.35100 40.08200
H  -4.56600 21.25200 38.85800
N  -10.36200 20.30600 42.60500
C  -10.01500 21.60400 42.35900
C  -11.22400 22.42800 42.22000
C  -12.33600 21.54600 42.38200
C  -11.71100 20.20500 42.59800
C  -11.19700 23.93200 41.97500
C  -13.79400 22.10700 42.18800
O  -13.91800 23.28100 41.86800
C  -15.04600 21.14100 42.31600
N  -10.60200 17.36400 42.77700
C  -11.94900 17.71100 42.85200
C  -12.88700 16.37600 42.87100
C  -11.75200 15.22000 42.76300
C  -10.48500 16.00900 42.84500
C  -13.77000 16.28100 44.15000
C  -11.80700 14.33600 41.44400
C  -11.81900 12.80500 41.57400
N  -7.76200 17.11800 43.06200
C  -7.94400 15.76000 43.00400
C  -6.61700 15.11000 42.96500
C  -5.67900 16.17500 43.03100
C  -6.45600 17.36700 43.03800
C  -6.31500 13.64100 42.85700
C  -4.23100 16.45300 42.98800
O  -3.23100 15.72200 43.03800
C  -4.16400 17.97000 43.08300
C  -3.54200 18.42400 44.28100
O  -4.15400 18.45500 45.38000
O  -2.21800 18.84600 44.04800
C  -1.49600 19.49100 45.14000
H  -8.63900 23.16900 42.16100
H  -13.51700 19.11600 42.86600
H  -9.44500 14.18400 42.98300
H  -4.61500 20.94700 43.53800
H  -5.92300 22.50900 41.35200
H  -5.58200 24.18000 43.08900
H  -5.32600 22.85000 44.27300
H  -6.97200 23.48100 43.91600
H  -3.14100 20.46400 41.91600
H  -3.60700 22.06200 41.29800
H  -5.03700 19.62900 40.18300
H  -3.36900 19.72900 39.82700
H  -11.41500 24.23200 40.95000
H  -10.31200 24.38600 42.42100
H  -11.94400 24.43600 42.58800
H  -14.90900 20.56900 43.23400
H  -15.05000 20.49700 41.43700
H  -16.00600 21.63900 42.18000
H  -13.59700 16.50100 42.05400
H  -11.83300 14.61300 43.66500
H  -13.23300 15.66500 44.87100
H  -14.78200 15.93400 43.94100
H  -13.98700 17.24300 44.61400
H  -10.98000 14.66700 40.81500
H  -12.76100 14.64800 41.01900
H  -11.51100 12.31300 40.65100
H  -12.86600 12.60200 41.79800
H  -11.12600 12.60300 42.39100
H  -7.19000 13.22900 42.35500
H  -6.33400 13.32500 43.90000
H  -5.34100 13.39500 42.43300
H  -3.56700 18.21200 42.20400
H  -1.27600 18.88800 46.02100
H  -2.03900 20.37600 45.47300
H  -0.61300 19.94900 44.69400
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


