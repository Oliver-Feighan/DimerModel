%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 35.66000 1.74600 29.23200
C  33.47200 2.34100 31.89500
C  38.12100 1.24100 31.62500
C  37.79000 1.61500 26.77900
C  33.01700 2.53400 27.01000
N  35.87600 2.18500 31.54100
C  34.82500 2.24700 32.36500
C  35.34000 2.18200 33.76400
C  36.88000 2.27300 33.63100
C  37.00600 1.84600 32.18200
C  37.59500 3.59000 33.93800
C  34.89600 0.87200 34.48900
C  34.26200 0.97200 35.91400
H  34.68800 2.14100 36.85600
N  37.69800 1.31100 29.19900
C  38.45500 0.95400 30.30700
C  39.74800 0.60100 29.81800
C  39.73600 0.82500 28.38300
C  38.39000 1.25600 28.03600
C  40.94400 0.20000 30.70900
C  40.89400 0.72300 27.34900
O  40.78100 0.94900 26.18600
C  42.22500 0.23800 27.88900
N  35.39600 2.02600 27.20700
C  36.47200 1.98000 26.39000
C  36.13800 2.22700 24.89300
C  34.55100 2.25100 24.96900
C  34.26500 2.24500 26.46600
C  36.85300 3.59900 24.36100
C  33.86400 1.14100 24.11300
C  33.41700 -0.12600 24.82000
N  33.60700 2.29400 29.34300
C  32.66300 2.48800 28.37100
C  31.34700 2.76700 29.04600
C  31.64200 2.71200 30.42700
C  33.04500 2.46100 30.58700
C  30.05700 3.20300 28.40600
C  31.07900 2.83900 31.76900
O  29.91300 3.09600 32.08200
C  32.20100 2.44500 32.72500
C  32.26400 3.54000 33.75100
O  32.77700 4.63100 33.58200
O  31.49900 3.13900 34.75300
C  31.03100 4.13600 35.70300
H  38.92600 1.22300 32.36200
H  38.44500 1.39800 25.93300
H  32.19900 2.76200 26.32300
H  34.93300 3.02400 34.32300
H  37.28500 1.57100 34.36000
H  36.96900 4.41100 34.28800
H  38.13600 3.93900 33.05800
H  38.32400 3.36100 34.71500
H  35.66700 0.10200 34.50700
H  34.14500 0.53000 33.77700
H  34.54900 0.03000 36.38300
H  33.17800 1.00900 35.80900
H  41.84000 0.79100 30.51900
H  41.29800 -0.80900 30.49700
H  40.73900 0.22800 31.78000
H  42.13800 -0.68600 28.46100
H  42.66900 1.05000 28.46500
H  42.94900 0.14100 27.08000
H  36.48000 1.41700 24.24900
H  34.11700 3.15500 24.54200
H  36.02900 4.30300 24.24400
H  37.25500 3.51600 23.35100
H  37.53000 4.05700 25.08300
H  34.65800 0.79900 23.44900
H  33.05000 1.52900 23.50100
H  34.25200 -0.36700 25.47900
H  33.30700 -0.80900 23.97800
H  32.56000 -0.11800 25.49300
H  30.04900 3.41600 27.33700
H  29.71100 4.11800 28.88700
H  29.33300 2.43700 28.68200
H  31.98900 1.57300 33.34400
H  31.22900 5.19300 35.52800
H  31.60100 3.89200 36.60000
H  30.07700 3.79700 36.10700
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


