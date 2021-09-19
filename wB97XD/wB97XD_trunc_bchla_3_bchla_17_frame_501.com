%nproc=24
%mem=175gb
#p wB97XD/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 1.29600 7.71300 26.88100
C  1.46900 9.82700 29.67000
C  2.42600 5.06900 28.94700
C  1.47300 5.72500 24.13400
C  0.92500 10.50000 24.87700
N  1.69900 7.49100 29.04900
C  1.51500 8.49600 30.00300
C  1.43800 7.84300 31.43000
C  1.83300 6.38800 31.15200
C  2.05200 6.30500 29.60100
C  3.12000 6.01700 31.93700
C  0.04000 8.06800 31.99300
C  -0.12000 8.13500 33.50300
H  1.05800 7.66300 34.38500
N  1.84800 5.68100 26.53500
C  2.31100 4.77600 27.50900
C  2.63500 3.49600 26.90400
C  2.31300 3.59100 25.49300
C  1.87500 5.04400 25.33500
C  3.03900 2.24400 27.69300
C  2.37500 2.55400 24.38800
O  2.03400 2.84400 23.24900
C  2.64700 1.14500 24.75200
N  1.13700 8.09300 24.78900
C  1.15400 7.08500 23.87300
C  0.77800 7.64200 22.47500
C  0.53200 9.18300 22.70900
C  0.76500 9.29400 24.23800
C  1.89300 7.33600 21.45900
C  -0.89600 9.67800 22.35700
C  -0.86300 10.38200 21.00400
N  1.30900 9.79300 27.10700
C  1.22300 10.82400 26.23200
C  1.50800 12.06900 26.98600
C  1.50300 11.75700 28.31400
C  1.42300 10.34100 28.32600
C  1.55300 13.42200 26.29600
C  1.55800 12.26300 29.69200
O  1.68500 13.41500 30.08900
C  1.44800 11.03500 30.57500
C  2.56000 11.13900 31.53300
O  3.76600 10.89700 31.38900
O  2.04800 11.65700 32.70500
C  2.99200 11.90900 33.77200
H  2.69600 4.26200 29.63100
H  1.45200 4.99400 23.32300
H  0.89300 11.33000 24.16800
H  2.21800 8.40600 31.94400
H  1.06000 5.77700 31.61800
H  4.00400 5.85500 31.32000
H  2.97100 4.99800 32.29300
H  3.38200 6.77100 32.67900
H  -0.56100 7.28300 31.53300
H  -0.37400 8.97400 31.55200
H  -1.02500 7.59600 33.78000
H  -0.31800 9.19800 33.64500
H  3.22100 2.42800 28.75200
H  3.95300 1.82100 27.27400
H  2.28900 1.47200 27.52500
H  1.81900 0.77800 25.35900
H  3.66100 1.02700 25.13400
H  2.61600 0.48900 23.88200
H  -0.10600 7.10300 22.13500
H  1.34500 9.68200 22.18000
H  2.33200 8.19000 20.94400
H  1.31000 6.77100 20.73200
H  2.77000 6.84600 21.88100
H  -1.17200 10.51700 22.99600
H  -1.64500 8.88600 22.34300
H  -1.71200 9.95300 20.47200
H  -0.02300 10.09700 20.37100
H  -0.94100 11.46900 21.03900
H  1.59500 13.38900 25.20700
H  2.38500 14.03800 26.63500
H  0.58900 13.91000 26.44400
H  0.47700 11.06900 31.06800
H  3.81000 12.49700 33.35500
H  3.40600 11.00100 34.21100
H  2.51900 12.68100 34.37800
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


