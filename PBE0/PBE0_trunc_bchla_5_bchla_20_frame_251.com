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


