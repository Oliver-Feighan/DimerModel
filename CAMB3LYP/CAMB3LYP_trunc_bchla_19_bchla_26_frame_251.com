%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 25.49300 50.21300 27.08700
C  23.37900 51.54000 29.65700
C  27.62800 49.17200 29.69500
C  27.85200 49.59100 24.90700
C  23.46200 51.62200 24.75200
N  25.44300 50.31400 29.39800
C  24.41900 50.80900 30.13400
C  24.66700 50.47600 31.59600
C  26.00400 49.68800 31.68900
C  26.38100 49.68400 30.19500
C  27.10100 50.37400 32.68300
C  23.49400 49.61900 32.23900
C  22.96400 50.09200 33.59300
H  23.92100 50.77900 34.64800
N  27.50700 49.52400 27.28400
C  28.12900 49.15000 28.38200
C  29.45200 48.73800 28.00600
C  29.55900 48.68700 26.59300
C  28.31600 49.30500 26.21200
C  30.44200 48.13300 29.00000
C  30.57100 48.09600 25.62400
O  30.49100 48.23300 24.41700
C  31.66500 47.26500 26.23600
N  25.76300 50.88600 25.17100
C  26.76300 50.29700 24.39500
C  26.51800 50.44800 22.86900
C  25.02700 50.95500 22.90400
C  24.70800 51.17300 24.38300
C  27.54000 51.40800 22.24300
C  23.97400 49.88700 22.34200
C  23.37400 50.12000 20.99600
N  23.75000 51.22900 27.14200
C  23.04900 51.77400 26.08800
C  21.90700 52.57500 26.52100
C  21.98600 52.51200 27.95800
C  23.17000 51.80300 28.27000
C  20.85200 53.08800 25.60900
C  21.31100 52.84500 29.20300
O  20.22900 53.44600 29.33300
C  22.22800 52.23700 30.41000
C  22.63900 53.39700 31.18500
O  23.23000 54.34500 30.71600
O  22.29300 53.27500 32.47200
C  22.54800 54.35700 33.47400
H  28.29000 48.76300 30.46100
H  28.36300 49.05600 24.10300
H  22.79700 51.98300 23.96500
H  24.78800 51.42300 32.12200
H  25.86900 48.70800 32.14700
H  26.61300 51.02100 33.41200
H  27.77100 50.92200 32.02000
H  27.67200 49.59700 33.19100
H  23.69200 48.55100 32.31900
H  22.62100 49.66400 31.58800
H  22.46400 49.28700 34.13200
H  22.16000 50.77100 33.30900
H  30.45600 47.04300 29.00900
H  30.32700 48.56400 29.99500
H  31.42700 48.48300 28.69200
H  31.35000 46.60400 27.04400
H  32.50000 47.83500 26.64100
H  31.97600 46.50800 25.51600
H  26.62800 49.45900 22.42500
H  24.97100 51.86600 22.30800
H  28.36600 51.68000 22.90000
H  27.08200 52.35600 21.96300
H  27.94500 50.95800 21.33600
H  23.09600 49.73600 22.97000
H  24.37200 48.87300 22.33000
H  23.66900 49.32700 20.31000
H  23.80100 51.07900 20.70200
H  22.28600 50.13300 21.06000
H  19.87400 53.08000 26.08900
H  20.94000 52.51200 24.68800
H  21.11300 54.12300 25.38400
H  21.65900 51.56400 31.05000
H  22.17700 55.30500 33.08400
H  23.60800 54.48100 33.69700
H  22.16500 54.13200 34.46900
Mg -9.26600 18.27700 42.96500
C  -5.81100 17.58000 42.98100
C  -8.47400 21.42900 42.13900
C  -12.47500 18.67000 42.48700
C  -9.75200 14.74100 43.14100
N  -7.34400 19.34600 42.56200
C  -6.06200 18.93000 42.69300
C  -5.04000 20.02200 42.45000
C  -5.98500 21.29000 42.45900
C  -7.37700 20.70000 42.36700
C  -5.82600 22.23200 43.72300
C  -4.08200 20.00700 41.18000
C  -4.36600 19.09100 39.96000
H  -3.96500 19.61900 38.53200
N  -10.31800 19.93500 42.50700
C  -9.83700 21.10100 42.09700
C  -10.97100 21.97200 41.68400
C  -12.19800 21.17200 41.95900
C  -11.67600 19.85500 42.32500
C  -10.82500 23.41900 41.29700
C  -13.61600 21.63300 41.74500
O  -13.81100 22.79100 41.38400
C  -14.89900 20.73300 41.88300
N  -10.94800 16.84400 42.71500
C  -12.12100 17.37600 42.70200
C  -13.18200 16.19300 42.82000
C  -12.25700 14.97500 42.77500
C  -10.85500 15.54600 42.86600
C  -14.16200 16.18100 44.00500
C  -12.30000 14.18000 41.36900
C  -12.88700 12.75700 41.40300
N  -8.10600 16.49700 42.98200
C  -8.41600 15.17700 43.17700
C  -7.19300 14.39200 43.48700
C  -6.09500 15.27800 43.29700
C  -6.75800 16.56600 43.10100
C  -7.09800 12.86900 43.72900
C  -4.68100 15.45100 43.26700
O  -3.71100 14.67400 43.28400
C  -4.42500 16.97400 43.13600
C  -3.61000 17.39300 44.36300
O  -2.45500 17.69200 44.51200
O  -4.50700 17.60300 45.41600
C  -3.95600 18.13500 46.70100
H  -8.18700 22.44600 41.86300
H  -13.54300 18.86100 42.35700
H  -9.93500 13.71100 43.45300
H  -4.44700 20.04900 43.36400
H  -5.75600 21.91400 41.59500
H  -5.11000 21.81600 44.43100
H  -6.79600 22.11100 44.20400
H  -5.58800 23.24500 43.39800
H  -3.08000 19.78600 41.54800
H  -3.99700 21.02300 40.79300
H  -5.44200 18.91700 39.97000
H  -3.77000 18.17900 39.94600
H  -11.77100 23.96000 41.28900
H  -10.25700 23.43000 40.36600
H  -10.11700 23.80900 42.02900
H  -15.08400 20.15000 42.78600
H  -14.85100 19.96100 41.11500
H  -15.81700 21.30800 41.76500
H  -13.87700 16.19700 41.98100
H  -12.39400 14.18200 43.51100
H  -14.02600 17.04700 44.65300
H  -14.09900 15.25100 44.57100
H  -15.20700 16.25900 43.70500
H  -11.30400 14.05300 40.94700
H  -12.88100 14.81200 40.69700
H  -12.44900 12.05100 40.69700
H  -13.92900 12.87300 41.10500
H  -12.71100 12.44000 42.43100
H  -7.98800 12.30200 43.45800
H  -6.94100 12.78400 44.80400
H  -6.18900 12.50600 43.25000
H  -3.84800 17.16600 42.23200
H  -2.93200 17.88100 46.97400
H  -4.58000 17.65400 47.45400
H  -4.05900 19.20000 46.91000


