%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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
Mg 6.73600 57.36100 41.87100
C  5.57600 54.03200 41.84400
C  9.54700 56.26000 40.29200
C  7.50300 60.62700 41.21800
C  3.48900 58.25400 43.03700
N  7.42000 55.35200 41.01200
C  6.89900 54.10900 41.23100
C  7.68800 53.00900 40.50000
C  9.08300 53.80000 40.39600
C  8.66300 55.23700 40.58100
C  10.16800 53.39100 41.42800
C  7.12500 52.63700 39.05500
C  7.43400 51.26500 38.47600
H  8.38100 51.07500 37.27200
N  8.34400 58.41300 40.90600
C  9.44600 57.74300 40.46000
C  10.47100 58.62900 40.12400
C  9.93400 59.94200 40.40300
C  8.52300 59.74000 40.83800
C  11.84700 58.22300 39.61000
C  10.74100 61.20700 40.21700
O  11.93000 61.13100 39.81300
C  10.22200 62.63600 40.60600
N  5.63100 59.10700 42.11600
C  6.19700 60.37100 41.79700
C  5.36400 61.47100 42.39400
C  3.99000 60.81500 42.65100
C  4.37000 59.30600 42.68800
C  5.80500 62.13900 43.71600
C  2.94100 61.09300 41.45000
C  2.15100 62.39300 41.46200
N  4.97100 56.31000 42.51500
C  3.75500 56.80200 42.96500
C  2.84300 55.78700 43.18500
C  3.53400 54.62400 42.80400
C  4.81800 55.00700 42.38100
C  1.45400 55.97300 43.70900
C  3.46600 53.25900 42.78000
O  2.61500 52.50900 43.21600
C  4.68000 52.84700 41.92100
C  4.09400 52.40300 40.61400
O  3.41000 53.09000 39.92100
O  4.40200 51.06200 40.37800
C  3.73000 50.54700 39.23600
H  10.54100 55.93300 39.97900
H  7.77500 61.68500 41.22200
H  2.57200 58.59300 43.52300
H  7.89200 52.14900 41.13700
H  9.49500 53.56500 39.41500
H  10.79300 52.70800 40.85200
H  9.75000 52.84400 42.27300
H  10.76800 54.24600 41.74000
H  7.52400 53.32000 38.30500
H  6.05400 52.83900 39.08200
H  6.45300 50.86600 38.21400
H  7.81100 50.65000 39.29300
H  12.72600 58.79500 39.90900
H  11.62200 58.24700 38.54400
H  12.11100 57.23800 39.99500
H  10.03600 62.67400 41.67900
H  9.32800 62.73100 39.99000
H  10.83200 63.47600 40.27400
H  5.23800 62.21700 41.61000
H  3.54200 61.14600 43.58800
H  5.15700 61.94700 44.57100
H  6.23700 63.12400 43.53800
H  6.70300 61.58800 43.99500
H  2.28600 60.23100 41.57500
H  3.43700 61.03400 40.48100
H  2.26900 62.89100 40.50000
H  2.53500 63.03900 42.25200
H  1.08300 62.20800 41.58000
H  0.80400 55.41700 43.03400
H  1.12000 57.01000 43.65700
H  1.32300 55.57400 44.71500
H  5.17600 52.01400 42.41800
H  2.78900 50.09000 39.54100
H  4.23500 49.69300 38.78300
H  3.52200 51.31600 38.49400


