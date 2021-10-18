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
Mg 17.10100 -2.07700 27.77000
C  16.45500 -0.19100 30.58100
C  19.23600 -3.94700 29.68100
C  17.96100 -3.74600 24.96300
C  15.37700 0.23000 25.82200
N  17.83400 -1.91300 29.91400
C  17.31200 -1.24300 30.86800
C  17.86400 -1.55200 32.35400
C  18.99500 -2.58100 31.95700
C  18.76000 -2.80300 30.38800
C  20.50900 -2.04800 32.15900
C  16.79200 -2.13400 33.37400
C  16.96200 -1.78000 34.89400
H  18.25700 -1.14000 35.29400
N  18.31900 -3.63600 27.38700
C  19.10500 -4.38300 28.35000
C  19.74900 -5.47800 27.63200
C  19.48000 -5.36700 26.24500
C  18.51900 -4.24500 26.14900
C  20.50800 -6.46400 28.44500
C  19.91900 -6.24400 25.04900
O  19.71600 -6.00400 23.89600
C  20.84100 -7.42500 25.41400
N  16.84100 -1.78700 25.63700
C  17.26600 -2.56700 24.69400
C  17.04700 -1.94200 23.28700
C  16.12300 -0.73100 23.56000
C  16.07300 -0.74300 25.10600
C  18.39300 -1.56900 22.63000
C  14.74200 -0.86000 22.84700
C  14.78100 -0.68700 21.32000
N  16.09100 -0.30700 28.09200
C  15.35600 0.45700 27.20300
C  14.73300 1.49900 27.97600
C  15.14100 1.35600 29.33600
C  15.94900 0.18700 29.29900
C  13.91400 2.52600 27.44500
C  15.02300 1.80700 30.68600
O  14.43000 2.83800 31.12300
C  15.98100 0.86900 31.51800
C  16.98700 1.71200 32.12800
O  17.99700 2.11700 31.55300
O  16.51500 2.07400 33.36400
C  17.31500 3.14200 34.00000
H  19.94300 -4.54500 30.25900
H  18.26400 -4.34800 24.10300
H  14.78900 0.88500 25.17500
H  18.31000 -0.60500 32.65900
H  18.90200 -3.48800 32.55400
H  20.47800 -0.95900 32.20400
H  21.12500 -2.36600 31.31800
H  21.03200 -2.40200 33.04800
H  16.90400 -3.21300 33.27100
H  15.74800 -2.07700 33.06400
H  16.83100 -2.63200 35.56100
H  16.22600 -1.02300 35.16500
H  21.52400 -6.11200 28.62600
H  20.42400 -7.48300 28.06700
H  20.00400 -6.52900 29.40900
H  20.33800 -8.26100 25.90100
H  21.56700 -6.93500 26.06200
H  21.42700 -7.84100 24.59500
H  16.51200 -2.73300 22.76200
H  16.63000 0.18800 23.26600
H  19.22100 -1.38900 23.31600
H  18.44300 -0.64500 22.05300
H  18.63700 -2.41000 21.98100
H  14.08400 -0.07400 23.21500
H  14.28900 -1.82100 23.09300
H  15.82700 -0.60000 21.02600
H  14.42000 0.31100 21.07100
H  14.22200 -1.38300 20.69500
H  13.25200 2.01200 26.74800
H  14.42200 3.30300 26.87400
H  13.44500 3.05700 28.27400
H  15.32300 0.48100 32.29600
H  17.17000 4.03200 33.38700
H  18.34600 2.83300 34.16900
H  16.87700 3.33200 34.98000


