%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 3.15800 -0.09400 44.31200
C  5.98100 1.96700 44.12700
C  1.59800 2.25100 42.30500
C  0.53500 -2.15100 44.13700
C  5.14900 -2.67000 45.58500
N  3.88700 1.89700 43.13200
C  5.03800 2.56500 43.33900
C  5.13100 3.86100 42.57000
C  3.59000 4.00400 42.27100
C  2.94600 2.64400 42.52000
C  2.95500 5.03300 43.24400
C  6.13600 3.75500 41.41600
C  6.21000 2.46100 40.57900
H  6.19500 2.64400 39.06600
N  1.34100 0.02700 43.38000
C  0.86600 1.12500 42.71300
C  -0.59800 1.03700 42.63500
C  -0.95100 -0.31800 43.02000
C  0.33900 -0.85500 43.55100
C  -1.47900 2.10800 41.88500
C  -2.32500 -0.98100 42.89500
O  -3.28500 -0.34000 42.56900
C  -2.49000 -2.43600 43.22700
N  2.98400 -2.23900 44.61900
C  1.73200 -2.88100 44.59300
C  1.80300 -4.39600 45.06100
C  3.31600 -4.49100 45.58100
C  3.86000 -3.01900 45.28400
C  0.72300 -4.97700 46.10400
C  4.23400 -5.64500 44.96900
C  5.28000 -6.17800 46.06800
N  5.16300 -0.35600 44.84000
C  5.79800 -1.43100 45.42800
C  7.17100 -1.05900 45.74700
C  7.33000 0.28900 45.33400
C  6.08400 0.65300 44.70700
C  8.16200 -1.87500 46.53500
C  8.23900 1.46700 45.15000
O  9.42700 1.57800 45.35600
C  7.35200 2.56500 44.42300
C  7.28300 3.74500 45.35300
O  8.21800 4.56900 45.51300
O  6.00100 3.85400 45.87200
C  5.72900 5.05600 46.58700
H  1.07500 3.03200 41.74900
H  -0.29600 -2.83600 44.31900
H  5.70000 -3.41500 46.16300
H  5.37600 4.70000 43.22100
H  3.29100 4.31200 41.27000
H  1.94000 4.68000 43.42800
H  3.04000 6.04200 42.84300
H  3.31500 5.08500 44.27200
H  7.12000 3.81500 41.88000
H  5.98200 4.60800 40.75500
H  5.33900 1.84600 40.80400
H  7.15200 1.96400 40.80900
H  -1.46500 2.03800 40.79700
H  -1.16500 3.10900 42.18100
H  -2.53000 2.14700 42.16900
H  -3.47600 -2.80900 42.94800
H  -2.33200 -2.65600 44.28300
H  -1.79700 -2.90200 42.52600
H  1.68700 -4.88400 44.09400
H  3.35000 -4.60600 46.66400
H  -0.10700 -5.37000 45.51800
H  0.37900 -4.20700 46.79400
H  1.20300 -5.75500 46.69800
H  4.88900 -5.25300 44.19100
H  3.62800 -6.44500 44.54300
H  6.30700 -6.00900 45.74400
H  5.13600 -7.24800 46.21500
H  5.31200 -5.69700 47.04500
H  8.31600 -2.82100 46.01700
H  7.75800 -2.08300 47.52600
H  9.13600 -1.39100 46.60500
H  7.97000 2.88400 43.58400
H  4.67500 5.07800 46.86300
H  6.08300 5.93000 46.04000
H  6.20600 5.06600 47.56700
Mg 40.74100 8.40300 29.69900
C  42.66000 10.02500 32.01600
C  38.73200 7.21500 32.07700
C  39.17900 6.75200 27.36300
C  43.28300 9.22800 27.34100
N  40.79300 8.45600 31.86800
C  41.50500 9.39200 32.54400
C  40.94800 9.63500 33.92000
C  39.85300 8.49700 34.00200
C  39.82700 7.94800 32.59200
C  40.10000 7.44100 35.10600
C  40.44100 11.18100 34.18500
C  40.87000 11.88200 35.48300
H  41.23900 10.96600 36.69200
N  39.13800 7.34400 29.69200
C  38.39900 6.91700 30.74900
C  37.34700 6.00300 30.37600
C  37.33700 6.02300 29.01400
C  38.65500 6.62800 28.64700
C  36.61000 5.29800 31.50000
C  36.33400 5.36600 28.04500
O  36.44200 5.38700 26.84100
C  35.20700 4.65100 28.73000
N  41.30000 7.91900 27.67200
C  40.36700 7.31400 26.91500
C  40.80200 7.26400 25.44700
C  41.95100 8.32700 25.38700
C  42.23200 8.49800 26.89400
C  41.21600 5.86100 24.82900
C  41.58000 9.72400 24.76100
C  40.48400 10.59900 25.31000
N  42.57200 9.48100 29.57000
C  43.55300 9.66900 28.59100
C  44.67900 10.46300 29.09900
C  44.36400 10.63700 30.45500
C  43.08700 10.02200 30.64800
C  45.85600 11.01300 28.33500
C  44.78800 11.17000 31.66800
O  45.80500 11.75000 32.01500
C  43.61600 11.01200 32.67100
C  44.29700 10.49000 33.87800
O  44.63700 9.31800 34.05700
O  44.61300 11.49900 34.74100
C  45.54500 11.24000 35.84300
H  37.99000 7.00700 32.85100
H  38.68100 6.27000 26.52000
H  43.95700 9.44200 26.50800
H  41.68000 9.37000 34.68300
H  38.87700 8.95900 34.14900
H  40.48500 6.53200 34.64400
H  39.14600 7.20600 35.57800
H  40.83000 7.87400 35.79000
H  39.35300 11.12200 34.15300
H  40.80700 11.84200 33.40000
H  40.04300 12.46600 35.88600
H  41.65400 12.60800 35.26800
H  37.22900 4.66200 32.13300
H  35.89000 4.51000 31.28000
H  36.12800 6.14400 31.99000
H  34.67500 5.32500 29.40100
H  35.72300 3.80100 29.17600
H  34.37800 4.27900 28.12900
H  39.88500 7.56300 24.93800
H  42.83200 7.86300 24.94300
H  41.11200 5.08600 25.58800
H  42.28100 5.83900 24.60100
H  40.70000 5.56000 23.91600
H  41.36300 9.40900 23.74000
H  42.43100 10.40400 24.73900
H  39.92000 10.04700 26.06300
H  39.77700 10.86400 24.52400
H  40.79900 11.50000 25.83600
H  46.78200 10.48700 28.56300
H  45.97200 12.07500 28.55200
H  45.60900 10.77100 27.30100
H  43.07400 11.94700 32.81400
H  46.05900 12.14500 36.16900
H  46.18500 10.38200 35.64300
H  45.04400 10.86200 36.73400


