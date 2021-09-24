%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg -2.04400 17.49100 26.73800
C  -2.15500 15.53600 29.75200
C  -2.48400 20.32400 28.82400
C  -2.10000 19.56800 23.93100
C  -2.02900 14.73500 24.85000
N  -2.50300 17.81300 29.09300
C  -2.37100 16.88400 30.06400
C  -2.66300 17.53000 31.46600
C  -2.96600 19.01100 31.04000
C  -2.55400 19.10600 29.55300
C  -4.45800 19.41500 31.30700
C  -1.49000 17.42200 32.44000
C  -1.82900 17.03100 33.84100
H  -0.68700 16.80800 34.84300
N  -2.05000 19.60600 26.40500
C  -2.24700 20.56500 27.39300
C  -2.18900 21.87200 26.76900
C  -1.98200 21.71100 25.40900
C  -2.11500 20.23400 25.17200
C  -2.46700 23.16900 27.47300
C  -1.85900 22.77200 24.24000
O  -1.68000 22.55800 23.06600
C  -1.78900 24.24000 24.60100
N  -2.21300 17.23900 24.73600
C  -2.21800 18.17300 23.73300
C  -2.50700 17.54400 22.35000
C  -2.23600 16.00100 22.60900
C  -2.12000 15.99700 24.18100
C  -3.92800 17.89500 21.80400
C  -0.96100 15.44100 21.90800
C  0.37100 15.96600 22.43900
N  -2.11100 15.55600 27.12300
C  -2.02800 14.45800 26.22200
C  -2.06000 13.19200 26.95100
C  -2.05400 13.57600 28.30400
C  -2.07200 14.99400 28.37800
C  -2.02800 11.85900 26.36100
C  -2.05700 13.00500 29.64700
O  -1.98500 11.83800 29.99400
C  -2.21400 14.25700 30.66700
C  -1.01500 13.99300 31.58100
O  0.19500 14.13800 31.33000
O  -1.45400 13.50000 32.72600
C  -0.50700 13.41700 33.91100
H  -2.54500 21.19300 29.48200
H  -2.11900 20.13500 22.99700
H  -1.94800 13.88000 24.17600
H  -3.47200 16.92500 31.87400
H  -2.32200 19.72000 31.56200
H  -4.61000 20.30200 31.92200
H  -5.09400 18.58600 31.61800
H  -5.01600 19.78600 30.44700
H  -1.29200 18.47500 32.64100
H  -0.61800 16.92800 32.01000
H  -2.49700 16.18800 34.02500
H  -2.37000 17.90200 34.21100
H  -2.54100 23.02800 28.55100
H  -3.37300 23.68600 27.15300
H  -1.69300 23.93100 27.38000
H  -0.93800 24.47000 25.24200
H  -2.76500 24.39500 25.06100
H  -1.83600 24.73500 23.63100
H  -1.82000 18.06600 21.68300
H  -3.08800 15.32400 22.54300
H  -3.89400 18.55600 20.93800
H  -4.67200 18.26800 22.50700
H  -4.29000 16.90500 21.52600
H  -1.14700 15.45600 20.83400
H  -0.92500 14.37100 22.11200
H  0.30800 16.76300 23.17900
H  0.88600 16.37900 21.57200
H  0.96700 15.25500 23.01000
H  -2.06800 11.03700 27.07600
H  -1.09600 11.82400 25.79700
H  -2.85300 11.79600 25.65000
H  -3.07800 14.15100 31.32300
H  -1.14400 13.30000 34.78800
H  0.04800 14.34000 34.07600
H  0.23900 12.62700 33.81700
Mg 25.94600 0.51000 29.41700
C  27.75500 0.04800 32.42200
C  23.08400 0.98300 31.28700
C  24.18800 0.65700 26.52400
C  28.83500 -0.18800 27.55800
N  25.56000 0.39300 31.59400
C  26.46800 0.41400 32.65000
C  25.73500 0.61600 34.03900
C  24.24300 0.62700 33.61200
C  24.25900 0.61600 32.06700
C  23.30400 -0.50900 34.16600
C  26.25400 1.80200 34.95100
C  26.76400 1.39400 36.30000
H  26.45300 2.38700 37.48900
N  23.95800 0.76900 29.03000
C  22.92100 1.01800 29.90600
C  21.67800 1.14300 29.16100
C  21.94400 1.05800 27.76500
C  23.43600 0.83800 27.70800
C  20.32000 1.49500 29.80100
C  21.00100 1.19500 26.57400
O  21.36800 1.18900 25.41700
C  19.44600 1.46900 26.69600
N  26.43900 0.16900 27.31600
C  25.54800 0.29300 26.34100
C  26.18900 0.06900 24.92500
C  27.73100 -0.01000 25.26000
C  27.69500 0.01500 26.81800
C  25.54500 -1.13400 24.17700
C  28.51900 1.10000 24.59300
C  29.85400 0.64900 24.00300
N  27.89800 -0.11200 29.79400
C  28.99000 -0.26300 28.96100
C  30.08800 -0.69100 29.77200
C  29.69800 -0.55100 31.12700
C  28.32100 -0.19300 31.11800
C  31.42300 -1.10700 29.28200
C  30.17600 -0.59300 32.52000
O  31.23000 -0.96300 33.02900
C  28.81000 -0.12800 33.44600
C  28.50200 -1.15100 34.39600
O  27.57900 -1.93000 34.29600
O  29.31100 -1.04800 35.50500
C  29.06500 -2.12600 36.49500
H  22.18200 1.26000 31.83700
H  23.66500 0.64600 25.56600
H  29.70300 -0.30600 26.90500
H  25.82900 -0.35200 34.53000
H  23.76600 1.56200 33.90600
H  22.35300 -0.03400 34.40900
H  23.77900 -0.95200 35.04100
H  23.10600 -1.31400 33.45800
H  25.42400 2.50500 35.02300
H  27.09400 2.25600 34.42500
H  27.85000 1.29700 36.26800
H  26.30600 0.41800 36.45800
H  20.08400 2.49100 29.42800
H  20.46000 1.63600 30.87300
H  19.49000 0.83300 29.55100
H  19.31200 2.38400 27.27400
H  18.87900 0.57700 26.96400
H  19.10200 1.53100 25.66400
H  26.08400 0.92600 24.26000
H  28.04800 -1.00900 24.96000
H  26.34500 -1.78700 23.82700
H  24.86300 -0.83400 23.38100
H  24.93000 -1.70700 24.86900
H  28.70300 1.93600 25.26700
H  27.98200 1.61200 23.79300
H  30.10400 -0.36600 24.31200
H  30.60900 1.32100 24.41200
H  29.90400 0.79600 22.92400
H  31.58400 -1.10700 28.20400
H  31.68500 -2.00300 29.84600
H  32.04200 -0.24000 29.51400
H  29.03700 0.80200 33.96800
H  29.24700 -3.09800 36.03600
H  28.05200 -1.96800 36.86500
H  29.74500 -2.08500 37.34500


