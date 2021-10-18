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
Mg -5.21000 24.54100 27.19500
C  -3.54100 26.38400 29.69300
C  -6.38000 22.49900 29.72600
C  -6.78300 22.59100 24.87600
C  -3.73200 26.39300 24.81700
N  -5.00900 24.47500 29.48500
C  -4.39300 25.40300 30.31600
C  -4.77800 25.15700 31.76100
C  -5.38000 23.70500 31.74000
C  -5.61200 23.56800 30.25200
C  -4.48700 22.58300 32.40100
C  -5.82800 26.14800 32.28500
C  -5.66000 26.57200 33.74800
H  -5.92800 25.52000 34.86000
N  -6.48000 22.85300 27.28100
C  -6.76500 22.19100 28.42600
C  -7.52200 21.02600 28.10100
C  -7.81700 21.04500 26.69400
C  -6.98500 22.14200 26.23200
C  -7.94700 20.03100 29.14500
C  -8.63000 20.10000 25.80900
O  -8.90700 20.38500 24.62300
C  -9.08700 18.72400 26.28000
N  -5.32000 24.55800 25.17000
C  -5.97300 23.61500 24.38200
C  -5.69900 23.91000 22.94300
C  -4.77200 25.17400 22.83700
C  -4.60800 25.45200 24.34000
C  -5.24900 22.67900 22.10000
C  -5.37300 26.38000 22.08300
C  -4.94400 26.45700 20.60500
N  -3.89500 26.01800 27.18100
C  -3.38100 26.77900 26.13200
C  -2.50100 27.90900 26.61000
C  -2.55000 27.71500 28.00000
C  -3.42600 26.62300 28.30100
C  -1.78700 28.93500 25.74900
C  -2.03300 28.27400 29.24000
O  -1.24900 29.20200 29.44300
C  -2.78100 27.47500 30.39000
C  -1.66500 26.98400 31.33300
O  -0.58100 26.54900 31.00300
O  -2.03300 27.27600 32.64500
C  -1.03300 26.85400 33.62600
H  -6.79600 21.92300 30.55500
H  -7.24600 22.04300 24.05300
H  -3.20100 26.86300 23.98700
H  -3.94200 25.15700 32.46000
H  -6.32400 23.59600 32.27400
H  -3.58700 22.93500 32.90600
H  -4.09000 21.81900 31.73300
H  -5.03500 22.07100 33.19200
H  -6.83700 25.75100 32.17600
H  -5.77900 27.08800 31.73500
H  -6.23400 27.46900 33.98000
H  -4.65000 26.93300 33.94100
H  -7.01100 19.85400 29.67500
H  -8.29300 19.05500 28.80500
H  -8.67800 20.40400 29.86200
H  -8.25700 18.16700 26.71500
H  -9.45700 18.28200 25.35500
H  -9.85300 18.80700 27.05100
H  -6.66900 24.06000 22.46800
H  -3.77200 24.98400 22.44500
H  -5.16100 21.82900 22.77700
H  -4.26700 22.88500 21.67500
H  -6.05700 22.54500 21.38000
H  -5.16500 27.26900 22.67800
H  -6.45400 26.27700 21.98600
H  -4.65800 25.46600 20.25200
H  -4.05700 27.05100 20.38400
H  -5.81000 26.77000 20.02100
H  -2.54700 29.61400 25.36200
H  -1.26100 28.53800 24.88100
H  -1.10500 29.45700 26.42000
H  -3.50500 28.04800 30.96900
H  -1.16000 25.83800 34.00100
H  -1.12400 27.59400 34.42100
H  0.01800 26.83100 33.33800


