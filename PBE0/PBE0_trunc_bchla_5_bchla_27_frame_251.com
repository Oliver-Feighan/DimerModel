%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 24.10500 -6.98700 45.76300
C  26.35700 -4.60000 44.67200
C  21.60100 -5.39700 43.81500
C  21.97300 -9.61100 46.32100
C  26.64500 -8.78600 47.26000
N  24.01000 -5.13100 44.52100
C  25.10000 -4.37600 44.14200
C  24.78600 -3.32200 43.01500
C  23.18500 -3.32200 43.04600
C  22.83700 -4.67700 43.88000
C  22.42500 -2.06700 43.62200
C  25.32600 -3.65200 41.66800
C  26.32400 -2.73600 40.98600
H  26.30600 -2.74400 39.44900
N  22.01900 -7.39700 45.22300
C  21.19700 -6.55000 44.45600
C  19.87600 -7.09000 44.39600
C  20.00000 -8.43900 45.00800
C  21.38600 -8.56500 45.55600
C  18.73000 -6.47100 43.58900
C  18.81400 -9.41900 45.03000
O  17.69600 -9.24400 44.55400
C  19.11500 -10.83500 45.37300
N  24.32600 -9.09200 46.54400
C  23.25000 -9.84400 46.81300
C  23.77300 -11.22400 47.20900
C  25.30200 -11.00000 47.62300
C  25.47100 -9.49100 47.13900
C  22.98600 -11.74600 48.46300
C  26.42600 -11.94400 47.02800
C  27.21200 -12.87800 47.93600
N  26.13500 -6.72100 46.02000
C  27.06100 -7.52500 46.69200
C  28.37300 -6.96800 46.56500
C  28.18100 -5.79800 45.76700
C  26.75600 -5.70300 45.52000
C  29.66500 -7.49100 47.10900
C  28.78000 -4.66700 45.18800
O  29.95100 -4.34500 45.21400
C  27.66100 -3.89000 44.39200
C  27.78800 -2.45400 44.68600
O  28.26300 -1.64500 43.89700
O  27.32400 -2.18500 45.96500
C  27.46400 -0.76100 46.34400
H  20.74200 -4.97500 43.29000
H  21.20500 -10.36000 46.52700
H  27.44000 -9.39000 47.70200
H  25.15500 -2.34400 43.32200
H  22.85800 -3.53200 42.02700
H  23.15600 -1.26300 43.70700
H  21.93600 -2.23300 44.58200
H  21.62700 -1.83100 42.91900
H  24.47200 -3.65600 40.98900
H  25.82900 -4.61800 41.69600
H  27.32300 -3.08700 41.24300
H  26.19000 -1.68500 41.24100
H  18.88600 -5.59200 42.96400
H  17.97200 -6.31900 44.35700
H  18.31100 -7.32700 43.06100
H  19.07300 -11.15900 46.41200
H  20.04000 -11.18300 44.91300
H  18.31700 -11.42600 44.92200
H  23.69500 -11.95400 46.40300
H  25.39200 -10.93900 48.70700
H  22.65700 -12.75800 48.22700
H  22.13300 -11.13900 48.76900
H  23.59300 -11.78000 49.36800
H  27.10900 -11.25400 46.53300
H  25.90100 -12.60500 46.33700
H  28.14800 -13.06200 47.40800
H  26.68200 -13.82800 48.00800
H  27.49600 -12.55500 48.93700
H  30.42900 -6.77800 46.79800
H  29.79900 -8.49600 46.71000
H  29.72100 -7.55200 48.19600
H  27.94500 -4.03400 43.34900
H  28.43100 -0.59300 46.81800
H  26.59100 -0.52900 46.95400
H  27.37100 -0.10100 45.48200
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


