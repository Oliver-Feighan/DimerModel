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
Mg 40.46600 8.46800 29.02600
C  42.13100 9.82400 31.82300
C  38.27200 6.93100 31.17900
C  38.87200 7.07400 26.41900
C  43.04200 9.56800 26.92900
N  40.22800 8.55200 31.19800
C  40.93100 9.11600 32.17600
C  40.40800 8.90400 33.61200
C  39.36400 7.73500 33.26600
C  39.31400 7.70200 31.75800
C  39.74100 6.39300 33.84500
C  39.69500 10.21300 34.26600
C  39.76000 10.23100 35.79600
H  40.52600 11.44400 36.38400
N  38.67500 7.29000 28.78400
C  37.94400 6.71700 29.85400
C  36.72400 6.18800 29.28700
C  36.85500 6.23500 27.90200
C  38.16400 6.90300 27.61900
C  35.64600 5.49400 30.09800
C  35.92100 5.81300 26.79100
O  36.21800 6.02700 25.65200
C  34.58600 5.19600 27.06600
N  40.97200 8.15700 27.00800
C  40.11100 7.61300 26.10000
C  40.78500 7.77300 24.70900
C  41.70800 8.97300 24.82300
C  41.98000 8.92500 26.34500
C  41.61300 6.52600 24.40200
C  41.12800 10.25100 24.34100
C  39.73300 10.63200 24.77100
N  42.22500 9.51500 29.25500
C  43.18900 9.89200 28.33900
C  44.24200 10.64200 29.04100
C  43.84000 10.65000 30.38800
C  42.65100 9.97700 30.46100
C  45.46200 11.19800 28.44400
C  44.09900 11.19900 31.69200
O  45.01500 11.97300 32.09800
C  42.93800 10.67900 32.66500
C  43.61600 9.87100 33.74100
O  44.43500 9.03300 33.50500
O  43.17500 10.12500 35.06800
C  43.60100 9.11400 36.08300
H  37.59700 6.39300 31.84700
H  38.44200 6.55400 25.56000
H  43.78300 10.08700 26.31700
H  41.23500 8.65200 34.27700
H  38.42900 8.08400 33.70600
H  39.16600 6.28900 34.76500
H  40.81800 6.26800 33.95900
H  39.27800 5.67800 33.16500
H  38.64200 10.24000 33.98400
H  40.12300 11.12800 33.85600
H  40.22100 9.26500 36.00300
H  38.77200 10.15700 36.25000
H  35.82300 5.83000 31.11900
H  35.93800 4.44600 30.03700
H  34.64100 5.77100 29.77900
H  34.06100 6.07300 27.44500
H  34.59100 4.42700 27.83800
H  34.14900 4.81100 26.14400
H  40.09900 7.86800 23.86800
H  42.63400 8.90500 24.25300
H  42.69600 6.55200 24.52300
H  41.40100 6.24300 23.37100
H  41.23400 5.71200 25.02100
H  41.17500 10.33500 23.25500
H  41.80700 11.00800 24.73300
H  39.85400 11.49300 25.42900
H  39.05700 9.95800 25.29800
H  39.19800 10.99600 23.89400
H  46.11400 10.36800 28.17200
H  45.85400 11.85300 29.22200
H  45.24700 11.72600 27.51400
H  42.43100 11.58700 32.99200
H  44.64200 9.27500 36.36300
H  43.48400 8.08500 35.74000
H  43.10400 9.17100 37.05100


