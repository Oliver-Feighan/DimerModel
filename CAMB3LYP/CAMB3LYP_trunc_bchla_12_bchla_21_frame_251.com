%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
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
Mg 15.77600 52.08200 25.40200
C  17.41000 50.69800 28.22000
C  13.20100 52.90800 27.42100
C  14.25400 53.22200 22.73700
C  18.25300 50.46100 23.37600
N  15.27500 51.74600 27.60200
C  16.22500 51.39100 28.54500
C  15.65700 51.73300 29.88000
C  14.15600 51.89200 29.63600
C  14.16000 52.23400 28.12500
C  13.39700 50.53700 29.89600
C  16.28100 53.03000 30.62900
C  16.58500 52.99000 32.15800
H  15.44800 52.47100 33.08500
N  13.93500 52.94500 25.15100
C  13.06900 53.20000 26.08700
C  11.87700 53.79900 25.46200
C  12.22000 54.11300 24.11200
C  13.51400 53.39600 23.93000
C  10.48000 54.11300 26.11100
C  11.50000 54.93700 23.04500
O  11.87100 55.02600 21.85900
C  10.19300 55.49000 23.36500
N  16.14700 51.69600 23.28800
C  15.39400 52.41700 22.36500
C  16.06600 52.28500 20.94400
C  17.38400 51.52800 21.24300
C  17.29200 51.23400 22.74200
C  15.13200 51.68100 19.89500
C  18.69000 52.21600 20.77200
C  19.80100 51.32100 20.17500
N  17.47400 50.82500 25.67600
C  18.34200 50.20900 24.78000
C  19.41200 49.59700 25.52600
C  19.05500 49.71600 26.87000
C  17.86300 50.48800 26.94600
C  20.58200 49.01900 25.00900
C  19.46500 49.44600 28.25600
O  20.44600 48.90400 28.74000
C  18.40300 50.10900 29.19300
C  17.85700 49.01900 30.09800
O  16.99100 48.18400 29.96600
O  18.53900 49.25300 31.25700
C  18.03900 48.60700 32.49400
H  12.40000 53.26000 28.07500
H  13.81800 53.79200 21.91400
H  19.12700 50.15700 22.79700
H  15.77300 50.97700 30.65700
H  13.74100 52.77300 30.12500
H  12.38500 50.66400 30.27900
H  13.84800 49.93600 30.68600
H  13.32800 50.06400 28.91700
H  15.69100 53.91000 30.37500
H  17.23100 53.29600 30.16400
H  16.68900 54.04900 32.39700
H  17.49400 52.38900 32.18400
H  10.24800 55.17100 26.23100
H  10.47100 53.63100 27.08900
H  9.61500 53.80800 25.52200
H  10.23700 56.09500 24.27100
H  9.45700 54.69600 23.49200
H  9.96700 56.22700 22.59500
H  16.38500 53.29100 20.67100
H  17.33800 50.51100 20.85400
H  15.24100 50.59600 19.90400
H  15.36700 51.95500 18.86700
H  14.10200 51.79300 20.23200
H  19.14000 52.87200 21.51700
H  18.48100 52.96500 20.00700
H  20.48500 51.15800 21.00800
H  20.30400 51.95600 19.44500
H  19.37300 50.39500 19.79200
H  21.13200 48.28300 25.59500
H  21.23200 49.84300 24.71300
H  20.26200 48.46700 24.12600
H  18.86400 50.90900 29.77200
H  18.09100 49.38900 33.25100
H  18.61600 47.75900 32.86200
H  16.99700 48.31500 32.36600


