%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 7.63600 56.80200 41.94500
C  6.62500 53.44600 41.38500
C  10.68100 55.97500 40.45800
C  8.03200 60.02500 41.40500
C  4.22800 57.38500 43.14900
N  8.45000 54.91300 40.78100
C  7.95700 53.66100 40.95000
C  8.99000 52.64800 40.49000
C  10.32900 53.33300 40.76400
C  9.81400 54.85400 40.60800
C  10.97600 53.14900 42.13800
C  8.68800 52.22300 38.98700
C  8.78600 50.76200 38.78500
H  9.47200 50.26400 37.49800
N  9.15000 57.91900 40.99500
C  10.37200 57.38400 40.51700
C  11.15700 58.50500 40.06500
C  10.49600 59.72000 40.38700
C  9.15900 59.29200 40.93300
C  12.51600 58.28800 39.40700
C  10.98000 61.08400 40.08900
O  12.05800 61.21700 39.52300
C  10.25300 62.36800 40.32000
N  6.19300 58.49600 42.17300
C  6.76600 59.71000 42.05200
C  5.81300 60.81300 42.54700
C  4.61500 59.97500 43.27000
C  5.02400 58.52600 42.89900
C  6.51300 61.95400 43.41400
C  3.16600 60.37300 42.84800
C  2.45300 59.49100 41.76700
N  5.85400 55.73100 42.22300
C  4.59600 56.10000 42.77900
C  3.77800 54.83100 43.01200
C  4.55200 53.84100 42.42800
C  5.79000 54.43400 41.99000
C  2.39800 54.58300 43.65100
C  4.62700 52.43800 42.17000
O  3.86200 51.52600 42.48800
C  5.87100 52.24300 41.37900
C  5.33600 51.93400 39.96200
O  4.91100 52.73800 39.13500
O  5.35200 50.52000 39.80700
C  4.79600 49.97200 38.63700
H  11.68700 55.76400 40.08900
H  8.33200 61.07300 41.34400
H  3.25200 57.43300 43.63500
H  9.01800 51.75000 41.10700
H  11.18700 53.14700 40.11900
H  10.60100 52.28700 42.68900
H  11.03000 54.09000 42.68600
H  12.00300 52.89400 41.87400
H  9.43500 52.68300 38.34000
H  7.77500 52.65300 38.57500
H  7.73900 50.46500 38.83400
H  9.33600 50.33800 39.62500
H  13.35000 58.53400 40.06400
H  12.53800 58.83400 38.46500
H  12.64500 57.23800 39.14400
H  10.18000 62.52800 41.39600
H  9.25600 62.43600 39.88300
H  10.86800 63.22900 40.06100
H  5.41300 61.22100 41.61900
H  4.75100 60.04000 44.34900
H  7.57600 62.02500 43.18400
H  6.34900 61.97700 44.49200
H  5.91600 62.83300 43.16600
H  3.09900 61.41000 42.52000
H  2.61700 60.24800 43.78100
H  1.86900 60.30400 41.33500
H  1.78000 58.70400 42.10800
H  3.21700 59.14200 41.07200
H  2.47400 53.89100 44.48900
H  1.98400 53.99700 42.83100
H  1.91500 55.55200 43.77500
H  6.40600 51.37000 41.75300
H  4.82800 48.88300 38.58500
H  5.36500 50.36200 37.79400
H  3.76000 50.31200 38.63400
Mg -10.08900 40.79200 42.94700
C  -8.61600 37.87300 41.73900
C  -7.78200 42.64100 41.17400
C  -11.69100 43.40600 44.11700
C  -13.01100 38.80800 43.84600
N  -8.29800 40.23100 41.61500
C  -7.89900 38.98100 41.23300
C  -6.67800 39.04800 40.38600
C  -6.44000 40.59600 40.17200
C  -7.57700 41.25700 41.04500
C  -4.98600 41.06900 40.42300
C  -6.83200 38.29900 39.04000
C  -5.67600 37.32300 38.82800
H  -4.88500 37.42500 37.48900
N  -9.75400 42.76600 42.70300
C  -8.72900 43.36600 41.97300
C  -8.87200 44.79700 42.20800
C  -10.00600 45.05500 43.11000
C  -10.52500 43.68700 43.42700
C  -8.10300 45.78800 41.48800
C  -10.49500 46.41500 43.39200
O  -9.97600 47.42800 42.97200
C  -11.74500 46.70800 44.14600
N  -12.01700 41.03600 43.77600
C  -12.30900 42.20100 44.43900
C  -13.68100 42.11000 45.07300
C  -14.24200 40.84200 44.54700
C  -13.08300 40.19800 43.96600
C  -13.70800 42.29100 46.64600
C  -15.26700 41.05700 43.38100
C  -16.32200 40.05900 43.24400
N  -10.72900 38.74300 42.86600
C  -11.88000 38.08300 43.31600
C  -11.64700 36.67300 43.29400
C  -10.42100 36.52600 42.65600
C  -9.90000 37.83800 42.36900
C  -12.60800 35.66600 43.80900
C  -9.50700 35.58200 42.14800
O  -9.55000 34.37200 42.15700
C  -8.30600 36.38400 41.59900
C  -7.07900 35.93500 42.27600
O  -6.65100 36.36600 43.36000
O  -6.55600 34.85100 41.60900
C  -5.28400 34.28000 42.02100
H  -6.98400 43.15900 40.63900
H  -12.35200 44.18100 44.51300
H  -13.89800 38.25100 44.15500
H  -5.86700 38.79100 41.06800
H  -6.60400 40.92400 39.14600
H  -4.42500 40.13600 40.35700
H  -4.81000 41.64900 41.32900
H  -4.66100 41.71400 39.60700
H  -6.76000 38.94200 38.16300
H  -7.79900 37.79800 39.08100
H  -6.24400 36.39600 38.75100
H  -5.06200 37.35000 39.72700
H  -7.39900 45.26300 40.84200
H  -7.46100 46.42800 42.09400
H  -8.77100 46.37100 40.85300
H  -11.59200 46.04300 44.99600
H  -12.59700 46.39100 43.54400
H  -11.78600 47.75300 44.45400
H  -14.19200 42.98500 44.67300
H  -14.60800 40.10900 45.26600
H  -12.74000 42.63900 47.00700
H  -14.11000 41.44200 47.20000
H  -14.30700 43.14200 46.96900
H  -14.72400 41.14600 42.44000
H  -15.78500 42.00900 43.49600
H  -17.15500 40.34300 43.88800
H  -15.96800 39.05000 43.45800
H  -16.80100 40.16400 42.27000
H  -13.09400 36.01500 44.72000
H  -12.17300 34.67900 43.96700
H  -13.39900 35.66300 43.05900
H  -8.06700 36.13500 40.56500
H  -5.23500 33.70700 42.94700
H  -4.52300 35.06000 42.03400
H  -5.07200 33.53500 41.25500


