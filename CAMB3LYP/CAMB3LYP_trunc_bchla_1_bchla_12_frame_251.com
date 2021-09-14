%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

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


