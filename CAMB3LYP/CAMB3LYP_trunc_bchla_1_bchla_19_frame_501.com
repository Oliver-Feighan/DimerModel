%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg -2.17200 16.97900 26.95200
C  -2.72400 15.15000 29.85500
C  -3.03700 19.79700 28.58300
C  -2.39300 18.59300 23.98800
C  -2.09700 13.93400 25.14200
N  -2.70100 17.33200 28.98300
C  -2.73500 16.51500 30.09200
C  -3.14200 17.32600 31.35300
C  -3.73400 18.61600 30.70600
C  -3.02700 18.64900 29.35900
C  -5.29100 18.59800 30.51100
C  -1.92300 17.38000 32.29200
C  -2.24700 17.17600 33.72300
H  -1.11900 16.77000 34.55900
N  -2.42900 18.95600 26.42700
C  -2.71700 19.97700 27.27500
C  -2.49100 21.22100 26.58300
C  -2.29300 20.85800 25.20200
C  -2.45300 19.40200 25.14600
C  -2.60700 22.55600 27.22700
C  -2.12700 21.73600 23.91900
O  -1.95600 21.25500 22.79700
C  -2.02300 23.23300 24.00700
N  -2.58000 16.29000 24.86100
C  -2.48200 17.22200 23.79200
C  -2.41100 16.57800 22.44700
C  -2.09800 15.11800 22.91000
C  -2.20000 15.11800 24.43300
C  -3.71000 16.74000 21.58700
C  -0.67500 14.61400 22.47000
C  0.50900 15.54800 22.82800
N  -2.25800 14.94800 27.37500
C  -2.16300 13.82100 26.54400
C  -2.25500 12.68900 27.38800
C  -2.44400 13.17700 28.67500
C  -2.44900 14.53100 28.60300
C  -2.11400 11.26600 26.95000
C  -2.70100 12.78200 29.98000
O  -2.81800 11.63600 30.45100
C  -2.89700 14.02900 30.89400
C  -1.76000 13.88300 31.95300
O  -0.54300 14.00600 31.72700
O  -2.34900 13.58900 33.17800
C  -1.46500 13.15200 34.30500
H  -3.42500 20.61700 29.19100
H  -2.30900 18.94200 22.95600
H  -1.87100 12.97600 24.66900
H  -3.88000 16.77500 31.93500
H  -3.40400 19.48400 31.27600
H  -5.81600 17.73100 30.91300
H  -5.54300 18.53800 29.45300
H  -5.71700 19.45700 31.02900
H  -1.55900 18.40500 32.23700
H  -1.19600 16.65200 31.93100
H  -2.95800 16.35000 33.74700
H  -2.76300 18.06300 34.09000
H  -2.53500 22.48300 28.31200
H  -3.46700 23.06800 26.79400
H  -1.74500 23.16300 26.95000
H  -2.35600 23.58400 23.03000
H  -0.98600 23.54300 24.13800
H  -2.68700 23.64600 24.76700
H  -1.57800 16.92500 21.83400
H  -2.85900 14.45200 22.50200
H  -4.54500 16.97100 22.24900
H  -4.01900 15.87700 20.99700
H  -3.69300 17.53600 20.84200
H  -0.75400 14.40000 21.40400
H  -0.59800 13.63700 22.94700
H  1.24700 15.06400 23.46800
H  0.07900 16.42200 23.31800
H  1.05200 15.89400 21.94900
H  -3.02600 10.76600 27.27600
H  -1.20700 10.86800 27.40500
H  -2.07500 11.10900 25.87200
H  -3.86000 13.93400 31.39500
H  -0.97000 12.24300 33.96400
H  -1.93600 12.93200 35.26300
H  -0.82600 14.03500 34.31800
Mg 25.47400 50.80800 26.36000
C  23.42800 51.77500 29.14000
C  27.78000 49.80000 28.46900
C  27.48600 50.74800 23.72700
C  22.89100 52.35600 24.24800
N  25.55500 50.75400 28.55700
C  24.64000 51.22200 29.52100
C  25.12200 51.01800 31.00700
C  26.56400 50.28500 30.68300
C  26.67300 50.24400 29.16200
C  27.79400 51.06900 31.31300
C  24.14500 50.22000 31.93900
C  24.75600 49.66700 33.22100
H  24.21600 50.30100 34.51600
N  27.35300 50.24000 26.12300
C  28.22400 49.80700 27.09500
C  29.50900 49.42000 26.55200
C  29.38500 49.70200 25.13400
C  27.97500 50.19600 24.94200
C  30.61900 48.81900 27.27800
C  30.38700 49.50900 24.00900
O  30.13300 49.66800 22.84800
C  31.78500 49.18000 24.46100
N  25.24600 51.52100 24.28400
C  26.25500 51.32600 23.39800
C  25.83300 51.76400 21.97000
C  24.32300 52.13900 22.16700
C  24.10500 51.95100 23.65600
C  26.64800 52.90800 21.27600
C  23.28700 51.29000 21.34900
C  23.20200 51.69600 19.86100
N  23.59000 51.86200 26.56100
C  22.66400 52.33500 25.66000
C  21.52300 52.77100 26.36200
C  21.75600 52.64400 27.68000
C  23.06700 52.09200 27.80400
C  20.26700 53.21900 25.62600
C  21.30400 52.86100 29.05800
O  20.28500 53.34300 29.48100
C  22.42600 52.30600 30.05900
C  23.02400 53.48800 30.75300
O  23.93000 54.23000 30.44700
O  22.33300 53.57900 31.94500
C  22.72800 54.66700 32.87100
H  28.62400 49.48900 29.08800
H  28.23800 50.84300 22.94000
H  22.18300 52.67600 23.48000
H  25.36700 52.00000 31.41200
H  26.54900 49.29900 31.14800
H  28.47500 50.35700 31.77900
H  27.60200 51.88000 32.01500
H  28.41300 51.49500 30.52400
H  23.78800 49.47000 31.23200
H  23.21100 50.69500 32.24000
H  25.83400 49.82400 33.19700
H  24.45000 48.62700 33.33100
H  30.63000 47.83900 26.80100
H  30.42700 48.80700 28.35100
H  31.61000 49.25400 27.15000
H  31.91600 48.09800 24.46700
H  32.15000 49.82600 25.26000
H  32.40300 49.70200 23.73100
H  25.75100 50.83500 21.40500
H  24.20500 53.21400 22.02900
H  26.08400 53.79400 20.98300
H  26.95200 52.49600 20.31400
H  27.50400 53.28000 21.83900
H  22.29300 51.43100 21.77200
H  23.59200 50.25400 21.50300
H  23.57200 52.71200 19.72200
H  22.15600 51.57400 19.58000
H  23.82700 51.09700 19.19800
H  19.37700 53.06800 26.23700
H  20.20500 52.68100 24.68000
H  20.36300 54.25800 25.31000
H  22.01700 51.57800 30.76000
H  21.79500 55.04700 33.28700
H  23.27800 55.48100 32.39900
H  23.31400 54.08500 33.58300


