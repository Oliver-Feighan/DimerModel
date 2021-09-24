%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 44.01800 2.91800 47.66100
C  42.06000 5.71500 47.01900
C  41.00600 1.01600 47.31200
C  45.79900 0.15700 47.76200
C  46.87900 4.91700 47.24100
N  41.83200 3.33000 47.09800
C  41.29400 4.59200 46.98100
C  39.79500 4.50000 46.68300
C  39.48700 3.02600 47.06600
C  40.84700 2.38200 47.12900
C  38.69200 2.74400 48.40100
C  39.44800 4.83900 45.22400
C  38.00300 4.93400 44.91400
H  37.59500 5.76800 43.64000
N  43.49300 0.83400 47.54100
C  42.21400 0.26700 47.42800
C  42.35900 -1.16200 47.43400
C  43.77600 -1.43800 47.42700
C  44.43200 -0.09700 47.60800
C  41.15600 -2.13400 47.28600
C  44.34200 -2.85900 47.34100
O  43.62700 -3.84400 47.42200
C  45.80100 -3.04100 47.16300
N  46.06200 2.62900 47.74100
C  46.53800 1.38200 47.82400
C  48.11500 1.40300 47.69600
C  48.51000 2.91300 47.47000
C  47.02000 3.56100 47.42600
C  48.93400 0.68700 48.80200
C  49.24600 3.07700 46.12300
C  50.25900 4.22700 46.06800
N  44.47400 4.92000 47.22100
C  45.65200 5.61300 47.20700
C  45.34700 7.04400 47.09000
C  43.98000 7.11700 47.10900
C  43.47900 5.83500 47.12000
C  46.36600 8.18800 47.10600
C  42.87400 8.02300 47.02800
O  42.81600 9.22700 46.97900
C  41.57700 7.14800 46.98100
C  40.79100 7.44500 48.21000
O  41.19100 7.46200 49.34100
O  39.43700 7.57100 47.85400
C  38.52600 7.78100 48.94300
H  40.02600 0.53600 47.28500
H  46.51600 -0.66500 47.82200
H  47.74600 5.57800 47.17700
H  39.30800 5.14400 47.41400
H  38.96600 2.43700 46.31200
H  37.88000 2.08800 48.08800
H  38.53600 3.57700 49.08700
H  39.33000 2.09900 49.00500
H  39.75400 4.10000 44.48300
H  39.93900 5.76800 44.93400
H  37.52300 5.35800 45.79600
H  37.59600 3.92700 44.82400
H  41.34100 -2.41600 46.24900
H  40.16900 -1.67200 47.25600
H  41.01800 -3.02000 47.90600
H  46.05900 -2.38200 46.33400
H  45.98200 -4.06600 46.83900
H  46.21200 -2.73700 48.12600
H  48.33800 0.83400 46.79300
H  49.07900 3.37000 48.27900
H  49.32800 -0.23900 48.38200
H  48.30000 0.56800 49.68000
H  49.82500 1.27000 49.03500
H  48.45300 3.22300 45.39000
H  49.73200 2.16800 45.76800
H  51.18600 3.72000 45.79800
H  50.32400 4.81900 46.98100
H  50.02100 4.90100 45.24500
H  46.22000 8.97700 46.36900
H  47.30400 7.65100 46.96300
H  46.44100 8.50600 48.14600
H  41.03900 7.44900 46.08200
H  38.18200 8.81200 48.86000
H  38.86300 7.60600 49.96500
H  37.69400 7.09100 48.80500
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


