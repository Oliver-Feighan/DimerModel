%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 29.10200 58.37000 41.70600
C  26.21100 56.85000 40.52500
C  31.00300 55.92400 39.98500
C  31.78200 60.36200 41.82300
C  27.09800 61.07100 42.76300
N  28.62400 56.74200 40.07100
C  27.36300 56.19900 40.01300
C  27.34200 54.94600 39.19500
C  28.87800 54.50900 39.24100
C  29.59900 55.79900 39.77000
C  29.21600 53.21300 40.13700
C  26.66600 55.21800 37.75400
C  26.06800 54.05700 37.13700
H  26.60500 53.76400 35.72000
N  31.13700 58.12700 41.18100
C  31.72600 57.03700 40.48100
C  33.15300 57.29600 40.39000
C  33.42500 58.57400 40.93600
C  32.09000 59.05300 41.42900
C  34.07200 56.38900 39.66000
C  34.85200 59.26600 40.92400
O  35.77500 58.62800 40.42000
C  35.15900 60.50100 41.76400
N  29.32700 60.43000 42.12600
C  30.57900 60.94600 42.21600
C  30.58100 62.56100 42.56700
C  29.07200 62.72600 43.08400
C  28.42300 61.36300 42.68500
C  31.60200 63.02200 43.67500
C  28.34100 63.90500 42.43900
C  27.75100 65.04300 43.30600
N  27.08200 58.74100 41.92000
C  26.43600 59.84700 42.39100
C  25.02600 59.57700 42.45300
C  24.91200 58.38300 41.64600
C  26.21500 57.96500 41.32600
C  23.95900 60.50800 42.95200
C  23.95200 57.50800 41.00600
O  22.77000 57.42900 41.15100
C  24.75000 56.37400 40.25200
C  24.55300 54.95400 40.82700
O  24.66200 54.66400 41.99900
O  24.13200 54.08300 39.80100
C  23.73600 52.76500 40.19600
H  31.53000 55.10300 39.49600
H  32.55200 61.12000 41.66600
H  26.51800 61.92400 43.12000
H  26.89400 54.10800 39.73000
H  29.20300 54.31600 38.21800
H  29.65600 52.49600 39.44400
H  28.33800 52.79800 40.63100
H  29.96300 53.53300 40.86400
H  27.50600 55.54800 37.14300
H  26.02700 56.07100 37.98400
H  24.99800 54.25700 37.09500
H  26.28800 53.16600 37.72500
H  33.60300 55.42900 39.44500
H  34.92200 56.13300 40.29300
H  34.37800 56.81900 38.70600
H  34.62000 61.39400 41.44600
H  36.23600 60.66600 41.71900
H  34.91600 60.25600 42.79800
H  30.70200 63.02600 41.58900
H  28.95600 62.73400 44.16800
H  32.36800 63.65500 43.22800
H  32.08100 62.14600 44.11400
H  31.02600 63.57600 44.41700
H  27.57000 63.56900 41.74600
H  28.92000 64.44500 41.69000
H  26.66500 65.11700 43.36900
H  28.07200 65.96200 42.81500
H  28.11700 64.98900 44.33200
H  23.07500 59.90100 43.14900
H  23.49000 61.12600 42.18700
H  24.26100 61.15900 43.77200
H  24.32100 56.40300 39.25000
H  23.95500 52.21700 39.27900
H  22.65500 52.67800 40.30400
H  24.24700 52.46200 41.10900
Mg -0.43700 44.18600 24.82500
C  1.43400 43.80100 27.72700
C  -3.18200 43.12800 26.59100
C  -2.25600 44.31000 21.89400
C  2.46600 44.70300 22.96600
N  -0.82100 43.63800 26.94100
C  0.04100 43.65300 27.97200
C  -0.60100 43.41700 29.35900
C  -2.06400 42.96100 28.94000
C  -2.04900 43.21000 27.35600
C  -2.28200 41.47000 29.38200
C  -0.53600 44.82000 30.15500
C  0.12400 44.79000 31.56500
H  0.16600 43.36300 32.22200
N  -2.40500 43.72200 24.25200
C  -3.41200 43.37700 25.17700
C  -4.69200 43.40700 24.50200
C  -4.45000 43.69100 23.15200
C  -2.97800 44.00900 23.04400
C  -5.95900 42.85900 25.17400
C  -5.47200 43.73000 21.98500
O  -5.19700 43.99600 20.80000
C  -6.95700 43.49800 22.31300
N  0.09900 44.39300 22.73100
C  -0.90800 44.46700 21.72100
C  -0.26100 44.51400 20.32100
C  1.23600 44.92600 20.70600
C  1.26600 44.70400 22.24100
C  -0.47400 43.18300 19.61900
C  1.52900 46.44600 20.35700
C  0.65700 47.66400 20.81600
N  1.53600 44.17200 25.19400
C  2.63200 44.38000 24.35100
C  3.83100 44.31600 25.12200
C  3.42700 44.13900 26.43300
C  2.02300 43.94900 26.43800
C  5.25700 44.36900 24.61300
C  3.89700 44.01900 27.75600
O  5.03000 43.96500 28.11500
C  2.62100 43.67600 28.63100
C  2.82000 42.26800 29.12000
O  2.22200 41.30700 28.68300
O  3.77300 42.30100 30.05900
C  4.37700 41.01200 30.29200
H  -4.07700 42.87100 27.16100
H  -2.79100 44.11000 20.96400
H  3.38600 44.85200 22.39600
H  -0.09900 42.67700 29.98200
H  -2.86900 43.57400 29.34600
H  -2.77400 40.83600 28.64400
H  -2.90900 41.64900 30.25500
H  -1.32000 41.10700 29.74300
H  -1.54200 45.23900 30.14700
H  0.04200 45.50300 29.53200
H  -0.52900 45.47000 32.11100
H  1.13900 45.18400 31.60700
H  -6.35300 42.07900 24.52200
H  -6.71900 43.64100 25.20200
H  -5.72100 42.49500 26.17300
H  -7.06800 42.43500 22.53100
H  -7.49100 43.66100 21.37700
H  -7.30300 44.26900 23.00100
H  -0.72000 45.37600 19.83600
H  1.94200 44.36800 20.09200
H  -1.38200 42.67100 19.93800
H  0.36900 42.51600 19.79800
H  -0.61300 43.32900 18.54800
H  1.65900 46.51600 19.27700
H  2.46900 46.62800 20.87800
H  -0.12600 47.82100 20.07400
H  1.34000 48.51400 20.83800
H  0.27600 47.62500 21.83700
H  5.47400 45.30000 24.08900
H  5.31100 43.53800 23.91000
H  6.03200 44.40300 25.38000
H  2.44000 44.27600 29.52300
H  5.42200 41.06200 30.59700
H  4.25200 40.18600 29.59200
H  3.93900 40.47400 31.13300


