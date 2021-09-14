%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
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
Mg 0.26000 43.34200 24.75600
C  2.25300 43.09900 27.61800
C  -2.44200 42.42000 26.69400
C  -1.57400 42.87700 21.91500
C  3.14600 43.60500 22.86000
N  -0.01100 42.92700 26.88700
C  0.88400 43.10000 27.86800
C  0.29200 43.09400 29.28500
C  -1.20300 42.68500 28.99500
C  -1.29800 42.75900 27.37000
C  -1.68200 41.32300 29.54200
C  0.44600 44.41100 30.04600
C  0.46000 44.23200 31.57500
H  0.30000 42.82500 32.16800
N  -1.74500 42.84700 24.34000
C  -2.74500 42.55800 25.29700
C  -3.99500 42.34000 24.59700
C  -3.74800 42.40500 23.17500
C  -2.33100 42.75300 23.08200
C  -5.28700 42.18000 25.35000
C  -4.75800 42.32200 21.94400
O  -4.44000 42.53400 20.79900
C  -6.18400 41.93600 22.27200
N  0.72800 43.16000 22.71000
C  -0.18900 43.01100 21.71100
C  0.47200 43.04900 20.28400
C  1.92800 43.34600 20.63900
C  1.92700 43.36400 22.15600
C  0.26000 41.87700 19.24300
C  2.38800 44.69400 20.03600
C  2.04300 46.00200 20.78200
N  2.32100 43.50900 25.09100
C  3.38200 43.61100 24.22900
C  4.61600 43.69300 24.93600
C  4.22700 43.48400 26.33500
C  2.89200 43.36300 26.33900
C  5.99200 43.88800 24.47400
C  4.71200 43.25500 27.68900
O  5.85600 43.34100 28.14800
C  3.47500 43.00100 28.60000
C  3.56600 41.71300 29.28000
O  3.33400 40.61600 28.78300
O  4.17100 41.89000 30.49900
C  4.75700 40.73100 31.20000
H  -3.23100 42.10200 27.37900
H  -2.18100 42.91300 21.00800
H  4.04900 43.78500 22.27300
H  0.67700 42.26400 29.87700
H  -1.92300 43.45600 29.26700
H  -2.57300 41.63600 30.08600
H  -0.94100 40.89700 30.21800
H  -1.87200 40.62800 28.72400
H  -0.32500 45.13200 29.77700
H  1.38300 44.89600 29.77100
H  -0.37400 44.87500 31.85700
H  1.32100 44.72000 32.03400
H  -5.00100 42.06900 26.39600
H  -5.81000 41.30000 24.97500
H  -5.96200 43.02800 25.22500
H  -6.37400 40.88800 22.50600
H  -6.81500 42.27400 21.45000
H  -6.38700 42.55900 23.14300
H  -0.03500 43.91100 19.85100
H  2.51500 42.52100 20.23700
H  1.17400 41.51200 18.77600
H  -0.34600 42.26100 18.42300
H  -0.30900 41.09200 19.74200
H  1.97300 44.71300 19.02900
H  3.46900 44.59400 19.93400
H  1.75200 46.69500 19.99300
H  2.89700 46.44600 21.29300
H  1.29400 45.88900 21.56500
H  6.47500 44.58100 25.16200
H  6.09900 44.26000 23.45500
H  6.57000 42.97800 24.63600
H  3.48200 43.84300 29.29200
H  4.58400 39.79600 30.66600
H  4.26600 40.64900 32.17000
H  5.80700 41.02000 31.22900


