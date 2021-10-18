%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 46.81700 25.04400 28.06300
C  47.43500 27.22400 30.70400
C  46.22200 22.54500 30.24300
C  46.68300 22.91500 25.44100
C  47.91200 27.51000 25.91600
N  46.85500 24.93200 30.27600
C  46.92700 25.99900 31.13500
C  46.69200 25.55600 32.52900
C  46.35500 24.06200 32.46000
C  46.45000 23.76800 30.90100
C  47.14000 23.08900 33.32600
C  45.54800 26.46400 33.25900
C  45.96500 27.28000 34.50500
H  45.04500 27.40000 35.69500
N  46.38000 23.04600 27.85500
C  46.17200 22.14900 28.89200
C  45.97200 20.84900 28.36300
C  46.10600 20.92300 26.95000
C  46.38200 22.34000 26.69200
C  45.54700 19.69500 29.35300
C  46.04800 19.73200 25.99500
O  46.24200 19.78300 24.77100
C  45.65900 18.41400 26.50400
N  47.46200 25.07000 25.97900
C  47.20100 24.13500 25.08900
C  47.75300 24.52400 23.69900
C  47.64000 26.09100 23.86200
C  47.64600 26.26000 25.37800
C  49.21100 24.11600 23.27000
C  46.37000 26.80300 23.19000
C  45.02600 26.38600 23.62800
N  47.60700 27.03800 28.18700
C  48.00800 27.95100 27.19200
C  48.42700 29.16900 27.80900
C  48.29400 28.93400 29.16300
C  47.68500 27.64800 29.33600
C  49.11900 30.26800 27.16400
C  48.30700 29.55700 30.40200
O  48.64500 30.69400 30.74700
C  47.85400 28.48300 31.49300
C  49.08400 28.17700 32.28200
O  50.16500 27.73400 31.95500
O  48.78800 28.53400 33.53300
C  49.86000 28.29900 34.55500
H  45.83500 21.89200 31.02800
H  46.47600 22.23000 24.61700
H  48.19300 28.30200 25.21900
H  47.64500 25.71200 33.03300
H  45.30400 23.92600 32.71400
H  47.87300 22.82100 32.56600
H  46.56600 22.26800 33.75600
H  47.69200 23.45800 34.19000
H  44.70400 25.82300 33.51200
H  45.17000 27.27500 32.63600
H  46.05600 28.31300 34.16800
H  46.96700 27.02300 34.84600
H  46.22400 18.85300 29.21000
H  44.59100 19.46800 28.88100
H  45.60100 19.88500 30.42500
H  46.32300 18.09200 27.30500
H  45.94900 17.70600 25.72800
H  44.60000 18.29800 26.73600
H  46.98600 24.30500 22.95600
H  48.53000 26.46300 23.35500
H  49.97100 24.18900 24.04700
H  49.40300 24.79200 22.43600
H  49.23500 23.08400 22.91900
H  46.32700 26.67600 22.10800
H  46.28900 27.84300 23.50600
H  44.88500 26.94500 24.55200
H  44.84400 25.33800 23.86800
H  44.17200 26.64500 23.00300
H  50.15800 29.93700 27.14300
H  49.06200 31.22600 27.68000
H  48.71000 30.52400 26.18600
H  47.00000 28.80500 32.08900
H  49.50700 27.43600 35.12000
H  49.99000 29.19400 35.16300
H  50.84800 28.01700 34.19000
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


