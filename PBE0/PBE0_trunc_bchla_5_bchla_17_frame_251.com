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


