%nproc=24
%mem=175gb
#p wB97XD/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 1.27700 7.97600 26.81300
C  1.70500 10.25900 29.27500
C  1.82700 5.41400 29.09300
C  0.85500 5.69300 24.34900
C  1.21800 10.51800 24.46400
N  1.59300 7.92000 28.89400
C  1.77500 8.96200 29.74400
C  1.81600 8.52800 31.18600
C  2.13100 6.98700 31.08800
C  1.79700 6.74800 29.57900
C  3.53500 6.42400 31.56400
C  0.41400 8.75400 31.91200
C  0.55800 8.89600 33.45600
H  1.71800 8.24400 34.22700
N  1.40400 5.81900 26.67300
C  1.63700 4.96400 27.78200
C  1.59900 3.60000 27.34800
C  1.32800 3.64600 25.93400
C  1.17100 5.11200 25.62400
C  1.72200 2.42000 28.24700
C  1.18400 2.49900 24.94000
O  1.05900 2.62000 23.68600
C  1.20700 1.06700 25.44700
N  1.13600 8.07200 24.71300
C  0.77100 6.99400 23.90500
C  0.28500 7.38700 22.48000
C  0.38900 9.00800 22.61900
C  0.90800 9.24400 23.99000
C  1.23400 6.77600 21.44100
C  -0.88200 9.73900 22.30100
C  -0.72400 10.77600 21.22600
N  1.52400 9.91600 26.73700
C  1.45400 10.87900 25.73800
C  1.65900 12.16800 26.28700
C  1.81400 12.01500 27.69100
C  1.62000 10.64100 27.94600
C  1.67300 13.48200 25.60500
C  2.15100 12.69400 28.92300
O  2.45500 13.83200 29.22300
C  2.00700 11.52600 30.01000
C  3.06400 11.59000 31.05300
O  4.23800 11.79300 30.87600
O  2.59100 11.43200 32.29900
C  3.63300 11.47700 33.36400
H  2.18200 4.70500 29.84400
H  0.85400 4.95200 23.54700
H  1.09900 11.38200 23.80700
H  2.67100 8.94100 31.72100
H  1.35600 6.44600 31.63000
H  4.07900 7.25400 32.01500
H  4.16900 5.94900 30.81700
H  3.32200 5.63000 32.28000
H  -0.18400 7.89400 31.61000
H  -0.01700 9.66700 31.50100
H  -0.39200 8.49700 33.81200
H  0.49500 9.96300 33.66800
H  0.81200 2.12400 28.76900
H  2.52200 2.67200 28.94200
H  2.18600 1.56000 27.76500
H  2.26000 0.80200 25.54100
H  0.75700 0.48100 24.64600
H  0.71900 0.92900 26.41200
H  -0.69500 6.92500 22.36000
H  1.13100 9.36700 21.90500
H  2.13800 6.31200 21.83600
H  1.41900 7.51100 20.65800
H  0.69100 5.93800 21.00300
H  -1.27500 10.30600 23.14500
H  -1.71400 9.09900 22.00700
H  -1.48700 10.70900 20.45100
H  0.24000 10.68600 20.72600
H  -0.68400 11.71900 21.77200
H  2.02800 13.48300 24.57400
H  2.51200 13.99500 26.07600
H  0.69500 13.95600 25.68900
H  1.08300 11.88900 30.46100
H  3.95200 12.49700 33.58300
H  4.44500 10.77600 33.17000
H  3.18800 11.10700 34.28800
Mg 48.00100 15.36700 27.57000
C  45.75200 15.79400 30.49800
C  50.04400 17.47200 29.24900
C  49.70300 15.15900 24.90400
C  45.19600 13.91000 25.93200
N  47.85300 16.51500 29.59600
C  46.96000 16.37600 30.64500
C  47.41800 17.01300 31.96200
C  48.55800 17.91800 31.40100
C  48.86700 17.29700 29.98400
C  48.22700 19.37400 31.30400
C  47.82300 16.04400 33.12300
C  47.60700 16.50100 34.52700
H  48.68000 16.07500 35.52400
N  49.72200 16.17900 27.12000
C  50.48900 17.00200 27.97100
C  51.76400 17.20600 27.34600
C  51.69700 16.38300 26.13900
C  50.37000 15.89100 26.04500
C  52.97000 17.98900 27.98500
C  52.83500 16.13800 25.19600
O  52.69500 15.47100 24.14200
C  54.21900 16.61300 25.45000
N  47.52500 14.80200 25.63500
C  48.43700 14.77900 24.66100
C  47.90800 14.19700 23.31100
C  46.65000 13.42700 23.86100
C  46.39500 14.11500 25.23000
C  47.50200 15.31900 22.30200
C  46.93500 11.88100 24.07400
C  47.49500 11.04500 22.90400
N  45.91900 14.88800 28.08600
C  44.96000 14.24100 27.31500
C  43.72700 14.13900 28.06000
C  44.00400 14.69600 29.31500
C  45.32400 15.22300 29.27200
C  42.47900 13.44100 27.63500
C  43.42900 14.97200 30.62500
O  42.35300 14.72000 31.15100
C  44.65500 15.63100 31.44600
C  44.12500 16.85400 32.02400
O  43.69700 17.84400 31.41000
O  43.76300 16.58100 33.31200
C  43.14800 17.60200 34.09800
H  50.71300 18.03200 29.90600
H  50.25500 15.06200 23.96700
H  44.37700 13.44100 25.38200
H  46.60300 17.64800 32.31000
H  49.48300 17.75400 31.95300
H  48.08000 19.59600 30.24700
H  49.08100 19.96700 31.63100
H  47.34400 19.72300 31.83900
H  48.89700 15.88200 33.03600
H  47.33700 15.08800 32.93000
H  46.66100 16.03500 34.80200
H  47.32400 17.55200 34.58700
H  53.79100 17.37200 28.34900
H  52.63200 18.67300 28.76300
H  53.29200 18.60300 27.14300
H  54.25300 17.68500 25.25400
H  54.84800 16.07400 24.74100
H  54.46700 16.31100 26.46700
H  48.71600 13.57900 22.91800
H  45.74000 13.52600 23.27000
H  48.21500 15.26400 21.47900
H  47.43200 16.32600 22.71200
H  46.52600 14.96500 21.96700
H  45.99300 11.38800 24.31300
H  47.54700 11.81000 24.97300
H  48.43400 10.64000 23.28100
H  47.70000 11.78200 22.12700
H  46.83000 10.23200 22.61300
H  42.78200 12.85500 26.76700
H  41.77900 14.22600 27.34900
H  42.00600 12.85500 28.42300
H  44.98300 14.94300 32.22500
H  43.06700 18.54100 33.55000
H  43.82900 17.76200 34.93400
H  42.21500 17.25400 34.53900


