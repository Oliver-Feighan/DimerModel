%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
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
Mg -5.36100 24.76000 27.02500
C  -3.64200 26.63800 29.43600
C  -6.40000 22.65900 29.54700
C  -6.94500 22.99100 24.73300
C  -4.07500 26.91800 24.56900
N  -5.09300 24.76200 29.22400
C  -4.33400 25.60400 30.03700
C  -4.64200 25.43000 31.57600
C  -5.37000 24.10100 31.51700
C  -5.66000 23.80300 29.98800
C  -4.43700 22.92100 32.05400
C  -5.65800 26.61300 31.90200
C  -5.45100 27.38800 33.27700
H  -4.92400 26.43700 34.32100
N  -6.35800 23.02700 27.11800
C  -6.68000 22.31300 28.19300
C  -7.48400 21.16300 27.87400
C  -7.81800 21.31800 26.51400
C  -6.97100 22.47200 26.02900
C  -7.92400 20.11100 28.82500
C  -8.85200 20.60400 25.66800
O  -9.21200 21.01500 24.55900
C  -9.41800 19.25800 26.12900
N  -5.43800 24.92500 24.91100
C  -6.26300 24.05500 24.23800
C  -6.35200 24.50400 22.78100
C  -5.50900 25.87000 22.72900
C  -5.01000 25.97700 24.14600
C  -5.90500 23.33400 21.83500
C  -6.23300 27.10700 22.22100
C  -5.47900 27.85000 21.07400
N  -4.10000 26.41500 26.99800
C  -3.63700 27.19900 25.89500
C  -2.75500 28.21500 26.34700
C  -2.74300 28.01700 27.73800
C  -3.53100 26.89000 28.05100
C  -2.03600 29.28900 25.59600
C  -2.22900 28.50400 28.97800
O  -1.46700 29.40000 29.25600
C  -2.89200 27.68300 30.12000
C  -1.83900 27.24500 31.07200
O  -0.96900 26.39800 30.89000
O  -2.04600 27.91300 32.21100
C  -1.18400 27.58800 33.32800
H  -6.90800 22.00000 30.25400
H  -7.48200 22.35700 24.02300
H  -3.70300 27.64900 23.84800
H  -3.68500 25.58900 32.07200
H  -6.31300 24.16400 32.05900
H  -3.44500 23.27200 32.33900
H  -4.33700 22.17800 31.26300
H  -4.99100 22.51200 32.89900
H  -6.69800 26.29000 31.93800
H  -5.59700 27.30300 31.06000
H  -6.41300 27.65600 33.71400
H  -4.98200 28.36000 33.12200
H  -9.00400 20.22200 28.93100
H  -7.47900 20.25200 29.80900
H  -7.70000 19.10000 28.48600
H  -9.79600 18.88600 25.17700
H  -10.18200 19.34000 26.90300
H  -8.70400 18.49400 26.43600
H  -7.35300 24.74000 22.42000
H  -4.63900 25.71000 22.09200
H  -6.65200 23.29700 21.04200
H  -5.78500 22.37600 22.34100
H  -5.04100 23.62000 21.23400
H  -6.46800 27.83800 22.99400
H  -7.10900 26.70500 21.71000
H  -4.60800 27.27900 20.75200
H  -5.00600 28.76500 21.43000
H  -6.13000 28.01900 20.21600
H  -2.85300 29.97800 25.38000
H  -1.55300 28.80400 24.74800
H  -1.32200 29.82700 26.21800
H  -3.47700 28.43300 30.65400
H  -0.13900 27.68100 33.03100
H  -1.23800 26.56800 33.71000
H  -1.41700 28.31200 34.10900


