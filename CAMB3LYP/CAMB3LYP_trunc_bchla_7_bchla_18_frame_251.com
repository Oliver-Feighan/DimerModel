%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
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
Mg 34.81300 49.43200 25.14500
C  35.01600 47.65700 28.22400
C  33.32400 52.02800 26.83800
C  34.65200 50.90000 22.18900
C  35.90600 46.46600 23.43400
N  34.14900 49.79800 27.34500
C  34.46100 48.93800 28.38200
C  33.98800 49.65100 29.69500
C  33.38500 51.08500 29.24000
C  33.61700 50.95700 27.73100
C  31.94200 51.39000 29.69300
C  35.01900 49.70600 30.89700
C  34.55900 49.07300 32.17800
H  33.11300 48.96700 32.58500
N  34.17200 51.35300 24.58700
C  33.58900 52.24400 25.43200
C  33.45200 53.46600 24.68200
C  33.72500 53.14400 23.33500
C  34.19500 51.72200 23.26900
C  32.90900 54.73200 25.36600
C  33.46900 53.99900 22.13100
O  33.64100 53.56300 20.96800
C  32.94500 55.37500 22.31400
N  35.18300 48.68800 23.08400
C  35.11600 49.51200 22.09700
C  35.57200 48.95700 20.64400
C  35.95000 47.52900 21.06600
C  35.68700 47.53300 22.55200
C  34.36900 49.06600 19.63700
C  37.40800 47.05400 20.73300
C  38.54700 47.67100 21.57900
N  35.21400 47.49800 25.67400
C  35.61800 46.42200 24.84900
C  35.92400 45.27900 25.72000
C  35.71000 45.72300 27.04700
C  35.31700 47.10200 26.97200
C  36.42700 43.88500 25.23900
C  35.78900 45.30100 28.34400
O  36.05400 44.22400 28.83800
C  35.40700 46.52600 29.22100
C  34.34200 46.13400 30.11500
O  33.13700 46.19100 29.84500
O  34.87800 45.46700 31.22300
C  33.86400 44.85200 32.11500
H  32.63300 52.79500 27.19500
H  34.43900 51.33900 21.21200
H  36.24500 45.54000 22.96500
H  33.26200 48.90600 30.01800
H  33.96500 51.87500 29.71500
H  31.90300 52.32700 30.25000
H  31.43700 50.62000 30.27600
H  31.39300 51.49300 28.75700
H  35.30500 50.73900 31.09600
H  35.87800 49.12400 30.56300
H  34.83800 49.76200 32.97600
H  35.03100 48.09400 32.25300
H  32.96100 54.66600 26.45300
H  31.87600 54.89800 25.06300
H  33.55400 55.53200 25.00400
H  33.66100 55.90400 22.94300
H  31.99400 55.37100 22.84800
H  32.75700 55.84600 21.35000
H  36.44300 49.49500 20.26900
H  35.19100 46.83200 20.71000
H  33.44700 48.97400 20.21100
H  34.41800 48.21400 18.95800
H  34.29100 49.96700 19.02800
H  37.57700 47.10800 19.65800
H  37.43700 45.97400 20.87400
H  38.40300 48.69400 21.92500
H  39.48100 47.54600 21.03000
H  38.64300 47.01100 22.44100
H  37.30600 43.53600 25.78100
H  36.62200 43.84300 24.16800
H  35.74400 43.11800 25.60500
H  36.33300 46.69700 29.76900
H  33.38300 44.00700 31.62300
H  33.04700 45.51600 32.39500
H  34.36200 44.67800 33.06900


