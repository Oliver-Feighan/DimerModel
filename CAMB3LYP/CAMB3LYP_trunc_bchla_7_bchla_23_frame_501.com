%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 25.26700 0.12000 29.48500
C  27.14100 -0.08700 32.36300
C  22.44200 0.65100 31.48400
C  23.33900 0.24400 26.71700
C  27.97000 -0.76800 27.63500
N  24.82700 -0.00300 31.72000
C  25.83000 0.13900 32.68700
C  25.26100 0.57400 34.05200
C  23.73000 0.68200 33.70700
C  23.61100 0.48500 32.21100
C  22.87000 -0.30200 34.49000
C  25.80700 1.87500 34.74000
C  26.64900 1.68700 36.11600
H  25.98700 2.43400 37.27800
N  23.21300 0.29800 29.14100
C  22.23900 0.64000 30.09300
C  20.97400 1.01400 29.45200
C  21.16000 0.82900 27.97700
C  22.63600 0.50000 27.85200
C  19.74900 1.44500 30.22800
C  20.23300 0.96600 26.86300
O  20.59800 0.66100 25.71800
C  18.79000 1.40700 27.03500
N  25.61300 -0.29200 27.55300
C  24.69700 -0.04800 26.51300
C  25.29200 -0.18000 25.05300
C  26.81800 -0.54700 25.43400
C  26.79800 -0.51500 26.97900
C  24.57500 -1.29000 24.13700
C  27.83300 0.46900 24.80600
C  29.13600 -0.06400 24.29800
N  27.18000 -0.44400 29.94400
C  28.21700 -0.73500 29.05300
C  29.45900 -0.85000 29.73800
C  29.08500 -0.49600 31.08700
C  27.69300 -0.34600 31.12700
C  30.82000 -1.11200 29.24500
C  29.58400 -0.41400 32.42400
O  30.70000 -0.50100 32.94800
C  28.28000 -0.28500 33.33500
C  28.18400 -1.48700 34.20500
O  27.52600 -2.44100 34.06000
O  28.92900 -1.30000 35.36700
C  28.48600 -2.07100 36.51000
H  21.58300 0.90400 32.10900
H  22.84400 0.32400 25.74700
H  28.87000 -0.90500 27.03200
H  25.38700 -0.24700 34.75700
H  23.31500 1.64400 34.00900
H  22.43700 -1.01900 33.79300
H  22.09100 0.33900 34.90100
H  23.35600 -0.79100 35.33400
H  24.95600 2.52500 34.94500
H  26.49900 2.40700 34.08700
H  27.58200 2.20000 35.88200
H  26.76100 0.62900 36.35200
H  19.34600 0.54600 30.69600
H  19.02800 1.84900 29.51900
H  19.97500 2.10600 31.06500
H  18.28300 1.10300 26.11900
H  18.67500 2.47500 27.22300
H  18.18400 0.86400 27.76000
H  25.32400 0.80000 24.57800
H  27.16700 -1.50900 25.05800
H  25.34100 -1.79600 23.54900
H  23.91500 -0.92500 23.35100
H  24.00900 -2.02000 24.71600
H  28.05800 1.32200 25.44600
H  27.45300 0.92700 23.89300
H  29.88300 0.35500 24.97200
H  29.54600 0.28700 23.35100
H  29.13600 -1.15300 24.35500
H  30.88400 -0.76300 28.21400
H  30.89900 -2.19900 29.25700
H  31.53400 -0.73600 29.97800
H  28.35300 0.59900 33.96800
H  28.85600 -3.07500 36.30500
H  27.39700 -2.05300 36.47000
H  28.92200 -1.70000 37.43800
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


