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


