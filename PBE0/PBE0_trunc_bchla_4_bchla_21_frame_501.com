%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 8.37300 3.02000 27.58300
C  9.66800 1.59900 30.62700
C  6.70500 5.26800 29.53700
C  6.93600 4.14100 24.76700
C  9.73100 0.40000 25.89700
N  8.26000 3.41400 29.92600
C  8.92300 2.70600 30.88600
C  8.72800 3.42100 32.16900
C  7.67000 4.51900 31.87000
C  7.43900 4.42400 30.33700
C  6.33800 4.27100 32.64700
C  9.97700 3.96400 32.81700
C  9.90300 4.49900 34.23800
H  10.91500 4.00800 35.18600
N  7.11500 4.54900 27.20100
C  6.43400 5.32900 28.17600
C  5.69100 6.36700 27.43700
C  5.70300 5.96800 26.10300
C  6.57300 4.78400 25.98300
C  4.94900 7.46900 28.16400
C  5.04400 6.62700 24.87200
O  4.92600 6.07300 23.76300
C  4.55800 8.00900 25.03700
N  8.26900 2.33700 25.58000
C  7.65300 2.93000 24.57800
C  7.79400 2.17100 23.23000
C  8.87800 1.15300 23.59100
C  9.00200 1.33000 25.08700
C  6.44500 1.50000 22.80300
C  10.20800 1.29100 22.88500
C  10.86600 -0.01000 22.39200
N  9.53800 1.41400 28.13000
C  10.02500 0.42800 27.30000
C  10.88100 -0.54000 28.05000
C  10.69600 -0.06300 29.33200
C  9.85900 1.05900 29.36100
C  11.78700 -1.59900 27.50900
C  11.07900 -0.33300 30.71000
O  11.83900 -1.20400 31.11800
C  10.43500 0.73200 31.62700
C  9.77600 0.01400 32.68800
O  8.68700 -0.60500 32.61700
O  10.41500 0.23000 33.87100
C  9.76200 -0.17400 35.09700
H  6.23400 6.04600 30.14100
H  6.34700 4.39700 23.88400
H  9.91800 -0.50400 25.31400
H  8.28800 2.65400 32.80500
H  8.01700 5.51600 32.14200
H  5.84500 5.03100 33.25400
H  6.57700 3.45100 33.32400
H  5.53900 3.81600 32.06100
H  10.35300 4.83500 32.28000
H  10.84000 3.32300 32.63900
H  8.94300 4.25600 34.69100
H  10.21000 5.54200 34.16500
H  3.93400 7.63800 27.80500
H  5.45900 8.43200 28.14900
H  4.76300 7.27700 29.22100
H  5.36300 8.60800 25.46300
H  3.75000 7.88400 25.75700
H  4.37900 8.39200 24.03300
H  8.13900 2.90800 22.50400
H  8.61200 0.10400 23.46400
H  6.50600 0.41300 22.84700
H  6.19800 1.73200 21.76700
H  5.70700 1.85800 23.52100
H  10.93600 1.82600 23.49400
H  10.03500 2.08200 22.15600
H  11.47800 0.25700 21.53100
H  10.15100 -0.81800 22.23400
H  11.45900 -0.28500 23.26400
H  12.83900 -1.31500 27.47600
H  11.33000 -1.91300 26.57100
H  11.63300 -2.52100 28.07100
H  11.28200 1.20800 32.12100
H  8.73800 -0.52700 34.97200
H  9.83700 0.61100 35.84900
H  10.27900 -1.00300 35.58100
Mg 15.62100 51.82400 26.04800
C  17.47400 50.60700 28.62800
C  13.11000 52.75700 28.20700
C  13.81500 52.82700 23.36600
C  17.99400 50.42600 23.83300
N  15.23600 51.49000 28.16800
C  16.26900 51.27600 29.02600
C  15.94600 51.61800 30.44100
C  14.40600 52.07300 30.26100
C  14.21100 52.11900 28.78100
C  13.38100 51.07100 30.96200
C  16.88200 52.82900 30.90600
C  17.11300 53.00000 32.40900
H  16.27700 52.29800 33.50800
N  13.77400 52.76600 25.77300
C  12.87200 53.04400 26.81300
C  11.70500 53.64500 26.21000
C  11.90800 53.69900 24.81600
C  13.23300 53.12400 24.62900
C  10.53200 54.15900 27.07700
C  10.86200 54.19000 23.68400
O  11.18100 54.33400 22.49700
C  9.41300 54.55600 23.96000
N  15.93700 51.69300 23.86100
C  14.99700 52.16200 23.00600
C  15.27300 51.71300 21.57200
C  16.82600 51.21000 21.68000
C  17.02200 51.18500 23.20400
C  14.21700 50.69700 20.93700
C  17.94900 52.21900 21.11400
C  19.01900 51.70400 20.16900
N  17.42400 50.77200 26.12100
C  18.25400 50.26800 25.19100
C  19.45700 49.67400 25.85400
C  19.17000 49.86100 27.22700
C  17.86900 50.41600 27.30000
C  20.69200 49.07400 25.28100
C  19.74000 49.61300 28.49800
O  20.83000 49.08800 28.80500
C  18.66300 50.08200 29.47300
C  18.25000 48.96100 30.28600
O  17.83600 47.88700 29.97900
O  18.30900 49.39800 31.60600
C  17.88900 48.52200 32.76600
H  12.38500 53.27700 28.83600
H  13.18500 53.06500 22.50600
H  18.65800 49.80900 23.22500
H  16.01600 50.87800 31.23800
H  14.23400 53.06400 30.68100
H  12.50100 51.52900 31.41300
H  13.92900 50.42000 31.64300
H  12.99700 50.42600 30.17100
H  16.51300 53.78500 30.53400
H  17.84200 52.64500 30.42300
H  17.07200 54.07500 32.58900
H  18.14900 52.74600 32.63000
H  10.69700 54.07400 28.15100
H  9.65000 53.59900 26.76800
H  10.36400 55.20000 26.80000
H  9.47200 55.48000 24.53400
H  8.95500 53.91200 24.71100
H  8.82100 54.48200 23.04800
H  15.27500 52.60400 20.94300
H  16.96000 50.24900 21.18300
H  14.69800 49.77000 20.62400
H  13.70900 51.16400 20.09400
H  13.48400 50.38900 21.68300
H  18.44000 52.77000 21.91600
H  17.40700 53.00700 20.59200
H  19.97700 51.93800 20.63200
H  18.94100 52.12500 19.16700
H  18.97400 50.61500 20.12300
H  20.76700 49.12700 24.19500
H  20.76600 48.01700 25.53800
H  21.55200 49.59200 25.70500
H  19.01800 50.87900 30.12700
H  17.50200 47.57000 32.40200
H  17.07000 48.99900 33.30600
H  18.78800 48.39000 33.36800


