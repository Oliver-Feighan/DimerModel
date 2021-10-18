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
Mg 29.46400 59.14500 41.03800
C  26.57800 57.38900 40.28900
C  31.34900 56.53100 39.78200
C  32.15000 61.18900 41.08600
C  27.35400 61.82000 41.94900
N  28.99200 57.29200 39.87100
C  27.77200 56.74800 39.89300
C  27.77900 55.38900 39.12700
C  29.32900 54.99700 39.21800
C  29.96700 56.34400 39.73200
C  29.58100 53.72700 40.06600
C  27.37400 55.59400 37.68500
C  26.46000 54.52700 37.10300
H  26.57800 54.02100 35.61600
N  31.53300 58.78800 40.65000
C  32.11000 57.71400 40.08400
C  33.51500 57.95800 40.03300
C  33.80600 59.24700 40.57400
C  32.47900 59.81600 40.84000
C  34.47700 56.89100 39.66100
C  35.16200 59.85800 40.59300
O  36.13000 59.21400 40.22400
C  35.39700 61.26400 41.13200
N  29.71300 61.28200 41.56500
C  30.94400 61.84600 41.41400
C  30.85900 63.33400 41.52000
C  29.29200 63.52200 41.79300
C  28.69900 62.13600 41.84900
C  31.79900 63.96200 42.56600
C  28.51800 64.48000 40.78300
C  27.64200 65.63900 41.31400
N  27.43700 59.54700 41.14000
C  26.74400 60.62700 41.58900
C  25.41200 60.28000 41.81900
C  25.25400 58.98400 41.34800
C  26.48900 58.61100 40.88000
C  24.32700 61.25500 42.36300
C  24.37200 57.92800 41.16800
O  23.15900 57.91800 41.34600
C  25.14900 56.75500 40.39700
C  25.25800 55.58200 41.18700
O  25.99300 55.31600 42.14100
O  24.33400 54.68200 40.66500
C  24.17500 53.35900 41.24400
H  31.95200 55.65500 39.53500
H  32.91900 61.96200 41.15300
H  26.80500 62.68200 42.33400
H  27.08300 54.71900 39.63100
H  29.81700 54.88800 38.25000
H  28.66900 53.37400 40.54700
H  30.34000 53.93800 40.82000
H  30.00600 52.97200 39.40500
H  28.16800 55.74800 36.95400
H  26.92700 56.58800 37.67000
H  25.44000 54.85100 37.30900
H  26.66400 53.69100 37.77200
H  34.81300 56.52500 40.63100
H  35.29700 57.28300 39.05900
H  34.02200 56.08300 39.08800
H  35.10600 61.93900 40.32700
H  36.47900 61.37300 41.20200
H  34.84400 61.42100 42.05800
H  31.14300 63.61000 40.50500
H  29.10900 64.04800 42.73000
H  32.59600 64.44300 41.99800
H  32.19100 63.23400 43.27700
H  31.26900 64.75100 43.09900
H  27.85100 63.86700 40.17700
H  29.18800 65.01800 40.11200
H  26.69000 65.11600 41.40400
H  27.54400 66.41700 40.55800
H  27.77900 66.00900 42.33000
H  24.06200 61.95200 41.56900
H  24.63400 61.90700 43.18100
H  23.42900 60.71200 42.65900
H  24.61600 56.61800 39.45600
H  23.13300 53.07400 41.10200
H  24.52900 53.25500 42.27000
H  24.72400 52.65900 40.61400


