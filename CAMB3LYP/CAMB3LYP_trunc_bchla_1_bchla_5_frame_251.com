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


