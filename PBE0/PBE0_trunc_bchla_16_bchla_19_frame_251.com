%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 40.61000 41.38800 26.95600
C  39.86800 43.56800 29.58900
C  41.56600 39.09100 29.09600
C  42.01100 39.61800 24.36600
C  40.07100 44.02500 24.67100
N  40.59100 41.34900 29.06300
C  40.20000 42.29100 29.98900
C  40.16900 41.71800 31.40000
C  40.93300 40.36900 31.20600
C  40.94200 40.20500 29.72600
C  42.27000 40.16100 31.88400
C  38.74600 41.61400 32.03600
C  38.70900 41.06400 33.46600
H  37.95300 41.90400 34.48500
N  41.75400 39.53700 26.76500
C  41.90400 38.67000 27.78800
C  42.49400 37.45600 27.29600
C  42.73100 37.64100 25.92300
C  42.16600 38.98800 25.63800
C  42.83400 36.26800 28.17100
C  43.30000 36.69700 24.91500
O  43.61700 36.97700 23.72800
C  43.50600 35.23000 25.27200
N  41.00200 41.78000 24.87100
C  41.64700 40.91600 24.05300
C  41.55500 41.44600 22.52200
C  40.51800 42.61700 22.64600
C  40.44500 42.80800 24.12700
C  42.99600 41.87100 22.03700
C  39.09100 42.15200 22.17700
C  38.49200 40.88200 22.76100
N  40.00500 43.44200 27.05400
C  39.84200 44.33300 26.02900
C  39.43100 45.53300 26.59500
C  39.42300 45.36800 28.01000
C  39.78400 44.06700 28.24300
C  39.15600 46.83600 25.93100
C  39.20700 45.98800 29.33700
O  38.93200 47.12500 29.70300
C  39.39600 44.78700 30.41700
C  40.24300 45.23600 31.55400
O  41.22600 46.02900 31.45500
O  39.66700 44.79200 32.66300
C  40.07000 45.50000 33.86300
H  41.79000 38.25100 29.75700
H  42.42000 39.05300 23.52500
H  39.79000 44.77400 23.92900
H  40.67200 42.47500 32.00200
H  40.36400 39.50500 31.55100
H  42.49900 40.80900 32.73000
H  43.06200 40.23200 31.13800
H  42.35700 39.09400 32.09000
H  38.09900 41.02000 31.39000
H  38.28900 42.60300 32.04400
H  39.66400 41.20800 33.97100
H  38.39800 40.01900 33.48100
H  42.71500 36.43800 29.24100
H  43.88400 35.97700 28.18900
H  42.26200 35.37200 27.92800
H  44.31400 35.17500 26.00100
H  43.97700 34.67500 24.46100
H  42.58300 34.72100 25.55200
H  41.23200 40.63200 21.87300
H  40.99300 43.46400 22.15100
H  43.79600 41.51800 22.68800
H  43.14700 42.94900 21.97800
H  43.27300 41.44900 21.07200
H  39.29200 42.04500 21.11100
H  38.30900 42.90300 22.28500
H  39.25300 40.32600 23.30900
H  37.94200 40.24300 22.07100
H  37.70100 41.27300 23.40100
H  38.42700 46.78300 25.12200
H  40.10600 47.12900 25.48500
H  38.86200 47.61200 26.63700
H  38.38700 44.62200 30.79400
H  39.19900 45.70100 34.48800
H  40.53400 46.47500 33.71600
H  40.74300 44.82600 34.39400
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


