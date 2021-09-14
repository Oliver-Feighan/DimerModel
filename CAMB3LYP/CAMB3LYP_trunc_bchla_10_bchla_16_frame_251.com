%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 40.54200 8.71200 29.57000
C  42.46900 10.35300 31.83700
C  38.24000 7.70200 31.93300
C  38.87500 6.96300 27.14400
C  43.10200 9.33000 27.08900
N  40.50700 8.83000 31.66000
C  41.28100 9.79100 32.39300
C  40.70700 10.01900 33.73400
C  39.62800 8.84800 33.81600
C  39.41100 8.39300 32.39000
C  40.09700 7.61800 34.74300
C  40.24300 11.48300 34.00600
C  40.64100 12.04200 35.33400
H  40.55300 11.05000 36.56000
N  38.75400 7.55700 29.53000
C  38.00600 7.28200 30.63600
C  36.89200 6.50300 30.21600
C  36.97400 6.32300 28.78000
C  38.23000 7.03300 28.36500
C  35.87900 5.90900 31.21900
C  36.08100 5.65900 27.77900
O  36.28700 5.46300 26.57600
C  34.68000 5.19900 28.33400
N  40.98100 8.10600 27.39000
C  40.09800 7.49500 26.64000
C  40.46000 7.54300 25.15700
C  41.81600 8.30500 25.16100
C  42.00400 8.59400 26.64900
C  40.56900 6.08800 24.47000
C  41.66000 9.63200 24.29500
C  40.74000 10.78700 24.76200
N  42.48900 9.50900 29.49800
C  43.34000 9.72700 28.43100
C  44.44800 10.55700 28.94600
C  44.20300 10.86700 30.25000
C  42.96700 10.18300 30.54000
C  45.69600 10.84800 28.18600
C  44.66700 11.56800 31.44700
O  45.66700 12.23600 31.65000
C  43.53700 11.23400 32.51700
C  44.17900 10.46300 33.60100
O  44.82400 9.40400 33.50900
O  43.90700 11.19500 34.73900
C  44.28800 10.53200 35.99300
H  37.46200 7.40700 32.64000
H  38.23300 6.48700 26.39900
H  43.84700 9.59300 26.33600
H  41.53500 9.82500 34.41700
H  38.70200 9.27400 34.20100
H  40.95100 7.85100 35.37900
H  40.49700 6.86600 34.06300
H  39.17400 7.30300 35.22800
H  39.15600 11.54000 33.94800
H  40.73700 12.04000 33.21000
H  39.99600 12.85500 35.66500
H  41.60800 12.52200 35.48500
H  35.14900 6.60800 31.62700
H  36.53200 5.46300 31.96900
H  35.20300 5.18600 30.76300
H  34.82100 4.32000 28.96300
H  33.98600 5.15200 27.49500
H  34.17100 5.96300 28.92200
H  39.71300 8.14400 24.63900
H  42.64100 7.73100 24.73900
H  39.67900 5.98600 23.84800
H  40.55600 5.31000 25.23400
H  41.39700 5.86900 23.79700
H  41.39900 9.47200 23.24900
H  42.67000 10.03800 24.36000
H  41.45000 11.59500 24.93900
H  40.14900 10.52600 25.63900
H  40.01100 11.10100 24.01500
H  45.97000 11.83700 28.55200
H  45.54800 10.83600 27.10600
H  46.48500 10.23100 28.61400
H  43.32700 12.24400 32.86800
H  44.09900 11.31400 36.72900
H  45.36500 10.37600 35.93100
H  43.84100 9.55300 36.16700
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


