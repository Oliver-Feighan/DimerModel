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
Mg 24.45900 -7.80400 45.73000
C  26.58600 -5.04500 44.57700
C  21.88700 -6.07400 44.19000
C  22.84100 -10.55200 45.75100
C  27.52200 -9.37000 46.76800
N  24.34300 -5.76600 44.35600
C  25.30500 -4.83000 44.04500
C  24.68800 -3.74900 43.08400
C  23.16700 -3.94300 43.21100
C  23.11900 -5.33100 43.95600
C  22.33200 -2.87300 44.03000
C  25.23700 -3.82200 41.61800
C  25.48800 -2.49900 40.92700
H  24.41700 -2.03400 39.89100
N  22.60400 -8.23500 45.19100
C  21.64700 -7.40100 44.68000
C  20.38100 -8.15500 44.46900
C  20.67900 -9.47200 44.90600
C  22.09500 -9.47000 45.26500
C  19.08200 -7.63000 43.83600
C  19.81800 -10.73300 44.78900
O  18.67600 -10.59500 44.35900
C  20.22700 -12.02800 45.16700
N  25.20700 -9.77600 45.95500
C  24.19400 -10.69000 46.16700
C  24.84300 -12.02000 46.71400
C  26.38100 -11.70100 46.72900
C  26.35300 -10.16200 46.56900
C  24.30300 -12.52000 48.08200
C  27.32200 -12.42400 45.68300
C  27.76100 -13.86800 45.98200
N  26.63200 -7.37200 45.76700
C  27.69000 -8.05400 46.33400
C  28.89800 -7.26500 46.33900
C  28.50100 -6.05100 45.66300
C  27.12300 -6.13300 45.25800
C  30.21600 -7.74700 46.86200
C  28.99300 -4.84100 45.17500
O  30.13700 -4.41700 45.08200
C  27.84800 -4.18900 44.34700
C  27.67100 -2.81900 44.95900
O  27.62100 -1.78000 44.29600
O  27.39300 -2.80300 46.38500
C  26.88100 -1.51100 46.83900
H  21.04400 -5.51700 43.77700
H  22.37000 -11.50900 45.98500
H  28.34500 -9.94300 47.20000
H  24.83100 -2.76600 43.53500
H  22.69900 -4.12200 42.24300
H  21.79800 -2.30300 43.27000
H  22.98800 -2.17400 44.54800
H  21.67600 -3.27700 44.80100
H  24.55300 -4.27200 40.89800
H  26.19500 -4.34200 41.59800
H  26.49400 -2.54600 40.51000
H  25.50300 -1.72600 41.69600
H  18.22200 -8.02700 44.37600
H  19.15100 -7.95200 42.79700
H  19.01200 -6.54700 43.92800
H  20.35400 -12.03800 46.25000
H  21.09800 -12.40800 44.63300
H  19.40000 -12.62300 44.77800
H  24.71500 -12.87500 46.05100
H  26.80800 -11.85300 47.72000
H  25.17900 -12.58000 48.72800
H  23.87600 -13.51200 47.93800
H  23.57100 -11.76000 48.35800
H  28.20300 -11.78200 45.65400
H  26.82800 -12.52100 44.71600
H  27.32600 -14.00800 46.97200
H  28.83900 -13.92900 46.13800
H  27.42200 -14.64300 45.29500
H  30.43000 -8.65700 46.30300
H  30.17700 -8.08000 47.89900
H  31.04400 -7.08900 46.59800
H  28.22500 -4.06500 43.33200
H  26.34600 -1.90700 47.70300
H  26.35400 -0.97600 46.04900
H  27.73900 -0.96400 47.23000


