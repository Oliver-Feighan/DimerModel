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
Mg 8.73000 48.00700 24.72600
C  6.54700 48.35400 27.54300
C  11.38500 48.83200 26.74200
C  10.78900 48.16000 22.05300
C  5.85600 48.14800 22.68600
N  8.91100 48.51400 26.88400
C  7.92400 48.53400 27.88200
C  8.57600 48.79800 29.21500
C  9.96600 49.43600 28.79200
C  10.09000 48.94100 27.32100
C  9.83300 50.99700 28.89900
C  8.62100 47.55100 30.16700
C  8.25400 47.92900 31.51300
H  7.91700 46.69900 32.49700
N  10.84000 48.27700 24.55100
C  11.74400 48.56500 25.49500
C  13.05300 48.54500 24.91800
C  12.92300 48.37200 23.55600
C  11.42000 48.19500 23.29600
C  14.33000 48.90700 25.74100
C  14.02400 48.24600 22.48300
O  13.70400 48.13100 21.32700
C  15.48100 48.17400 22.87500
N  8.33400 48.41200 22.66100
C  9.38600 48.28200 21.76100
C  8.87500 48.15200 20.36300
C  7.28100 48.12500 20.55900
C  7.12600 48.16700 22.06100
C  9.40100 49.20400 19.44600
C  6.46300 46.89200 19.88400
C  5.44500 47.20900 18.86100
N  6.61400 48.15800 24.98900
C  5.61400 48.20500 24.04400
C  4.27700 48.23000 24.78200
C  4.59000 48.17400 26.15500
C  6.02700 48.18500 26.22900
C  2.90700 48.32700 24.09700
C  3.98700 48.26800 27.45000
O  2.78100 48.32100 27.78100
C  5.21700 48.56000 28.40400
C  5.00500 47.63900 29.64300
O  5.64000 46.57500 29.76100
O  4.15000 48.22900 30.51400
C  3.80400 47.24100 31.59300
H  12.25300 48.87100 27.40400
H  11.43600 48.18800 21.17300
H  4.90900 48.09600 22.14500
H  8.00400 49.60800 29.66700
H  10.77600 49.10700 29.44200
H  9.23900 51.44000 28.10000
H  10.83100 51.43100 28.94300
H  9.35800 51.24300 29.84900
H  9.67300 47.28200 30.26500
H  8.05700 46.70900 29.76500
H  7.29600 48.44900 31.52800
H  9.07600 48.51700 31.92200
H  14.96900 48.03700 25.59000
H  14.16900 49.02000 26.81300
H  14.74800 49.84700 25.38100
H  16.11600 48.17200 21.98900
H  15.42200 47.21200 23.38400
H  15.69700 49.01800 23.53000
H  9.23500 47.16500 20.07200
H  6.89900 49.04100 20.10900
H  10.07700 48.66800 18.78100
H  10.03800 49.93000 19.95100
H  8.64800 49.76300 18.89000
H  5.91200 46.49600 20.73700
H  7.19200 46.14700 19.56800
H  4.44900 46.91600 19.19200
H  5.77800 46.47600 18.12700
H  5.55100 48.21100 18.44300
H  2.66500 49.38900 24.07600
H  2.13000 47.78700 24.63800
H  2.97200 47.86900 23.11000
H  5.09400 49.61500 28.64500
H  3.12300 46.44700 31.28700
H  3.39400 47.68200 32.50100
H  4.75700 46.82700 31.92300


