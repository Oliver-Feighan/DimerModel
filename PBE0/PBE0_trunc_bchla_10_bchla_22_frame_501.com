%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg 9.20000 47.36600 24.88000
C  7.02600 47.76900 27.74700
C  11.57000 48.44700 26.79100
C  10.94600 47.19300 22.09000
C  6.29100 46.96700 22.73800
N  9.27600 47.92100 27.10800
C  8.32000 48.03400 28.09600
C  8.93500 48.52700 29.35600
C  10.25500 49.18700 28.88100
C  10.38000 48.52600 27.48800
C  10.19000 50.72500 28.75600
C  9.27100 47.32000 30.28400
C  8.91000 47.48600 31.81100
H  8.76300 46.21000 32.73200
N  11.12300 47.75000 24.52700
C  11.99300 48.08000 25.44400
C  13.30700 48.11300 24.83000
C  13.09300 47.62600 23.51900
C  11.67000 47.53200 23.33600
C  14.56100 48.44000 25.66100
C  14.15200 47.33300 22.49300
O  13.84400 46.95300 21.40300
C  15.59100 47.42400 22.88400
N  8.59400 47.38800 22.62900
C  9.59400 47.15300 21.78200
C  9.10700 46.72400 20.44900
C  7.56400 46.85800 20.64800
C  7.41500 47.09900 22.12700
C  9.70100 47.44000 19.18300
C  6.80800 45.57800 20.11700
C  6.00900 45.65200 18.78600
N  7.17400 47.35700 25.08600
C  6.10400 47.11400 24.19600
C  4.83100 47.11900 24.94600
C  5.17200 47.34400 26.27300
C  6.54300 47.52700 26.33600
C  3.42100 47.05000 24.31100
C  4.55900 47.48400 27.63300
O  3.43200 47.38700 28.05800
C  5.73500 47.94200 28.63300
C  5.68600 47.23700 29.86100
O  6.32700 46.19300 30.04600
O  4.84600 47.74900 30.84500
C  4.90400 47.16700 32.22100
H  12.46500 48.75400 27.33500
H  11.53100 47.07200 21.17600
H  5.40800 46.83900 22.11000
H  8.32300 49.32300 29.78100
H  11.13800 49.00400 29.49300
H  10.92100 51.19300 29.41600
H  9.18200 51.02300 29.04300
H  10.33700 51.10600 27.74500
H  10.35100 47.19900 30.19800
H  8.78700 46.45500 29.83100
H  8.00100 48.08400 31.87000
H  9.74500 48.08400 32.17800
H  15.18900 47.54900 25.70000
H  14.38200 48.83100 26.66200
H  15.20000 49.13500 25.11700
H  15.85600 46.71400 23.66700
H  15.91500 48.35300 23.35300
H  16.20400 47.15000 22.02500
H  9.34100 45.66500 20.34000
H  7.22100 47.78300 20.18400
H  10.49600 48.09000 19.54700
H  9.00000 47.94700 18.51900
H  10.12400 46.67900 18.52700
H  6.07000 45.27200 20.85900
H  7.53000 44.76900 20.22800
H  6.23700 44.78300 18.16900
H  6.25500 46.57800 18.26600
H  4.93600 45.62800 18.98000
H  3.49200 46.43400 23.41500
H  3.14200 48.04200 23.95500
H  2.68000 46.64400 25.00000
H  5.61100 49.00700 28.82700
H  5.88500 47.32300 32.67200
H  4.47300 46.16700 32.25100
H  4.32200 47.83100 32.86000


