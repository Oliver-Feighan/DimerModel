%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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
Mg 29.51000 58.90800 41.21400
C  26.38900 57.48600 40.17300
C  31.07800 56.05600 39.82600
C  32.44700 60.51200 41.46200
C  27.84100 61.85400 41.88400
N  28.82700 56.95900 40.12100
C  27.51800 56.70000 39.82500
C  27.43600 55.53200 38.85500
C  28.91100 54.97600 39.03000
C  29.71500 56.02800 39.69400
C  29.08100 53.62400 39.80100
C  26.96300 56.03100 37.45200
C  25.85600 55.22200 36.65900
H  26.19700 54.60200 35.29200
N  31.43900 58.34900 40.91600
C  31.91100 57.13400 40.35300
C  33.36200 57.11100 40.33800
C  33.74700 58.43000 40.81400
C  32.51100 59.13800 41.09000
C  34.11400 55.98100 39.76900
C  35.21500 58.86100 40.97700
O  36.08800 58.03800 40.57400
C  35.69500 60.19600 41.58500
N  30.14600 61.00100 41.36300
C  31.38600 61.37000 41.62100
C  31.40400 62.84700 41.84900
C  29.90900 63.22000 42.03600
C  29.21400 61.95000 41.65200
C  32.31600 63.46200 42.99500
C  29.40100 64.48800 41.17400
C  28.93400 65.75600 41.93000
N  27.56300 59.51200 41.26600
C  27.03200 60.71300 41.55100
C  25.61900 60.69600 41.44600
C  25.33000 59.42300 40.94900
C  26.56000 58.73200 40.80800
C  24.61100 61.72700 41.76100
C  24.21900 58.63400 40.51400
O  23.05400 58.95300 40.43300
C  24.92200 57.26000 39.91000
C  24.42400 56.08500 40.63500
O  24.27500 56.00600 41.86100
O  24.19400 55.10700 39.69100
C  23.68300 53.83100 40.18900
H  31.62000 55.20100 39.41900
H  33.35900 61.11200 41.46000
H  27.25300 62.67600 42.29800
H  26.78700 54.79200 39.32400
H  29.19500 54.98000 37.97700
H  29.75800 53.81000 40.63600
H  29.53300 52.82700 39.21100
H  28.17800 53.23200 40.27000
H  27.91200 55.97700 36.91900
H  26.65200 57.07400 37.39700
H  25.02500 55.92000 36.55900
H  25.42300 54.48900 37.33900
H  33.71100 55.81900 38.76900
H  33.89700 55.07100 40.32800
H  35.20300 56.02600 39.79800
H  35.47200 60.98900 40.87100
H  36.68900 59.96900 41.96900
H  35.05300 60.30300 42.46000
H  31.80800 63.29500 40.94200
H  29.82800 63.26800 43.12200
H  32.55800 64.46300 42.63700
H  33.18900 62.82200 43.12200
H  31.69200 63.52000 43.88700
H  28.59500 64.07600 40.56700
H  30.17400 64.71200 40.43900
H  27.96900 66.12300 41.58300
H  29.69100 66.53100 41.80700
H  28.89100 65.63000 43.01200
H  23.61700 61.29700 41.88000
H  24.71500 62.31200 40.84700
H  25.04500 62.40800 42.49300
H  24.83500 57.21500 38.82500
H  22.80800 53.96700 40.82500
H  24.47100 53.30400 40.72900
H  23.37200 53.22500 39.33800


