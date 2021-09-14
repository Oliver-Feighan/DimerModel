%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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


