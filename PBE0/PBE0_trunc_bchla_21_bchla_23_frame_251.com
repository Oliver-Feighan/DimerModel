%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
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
Mg -9.42700 40.51200 41.96000
C  -8.27600 37.27100 41.01300
C  -7.18900 41.70700 39.82100
C  -11.06300 43.46900 41.99400
C  -12.24300 38.82400 43.34600
N  -7.88600 39.63200 40.71800
C  -7.50800 38.28500 40.52900
C  -6.30600 38.11700 39.54700
C  -5.96100 39.60100 39.23600
C  -7.01600 40.39000 40.04000
C  -4.48700 39.82100 39.61500
C  -6.54600 37.26600 38.22600
C  -5.77300 35.95800 38.10300
H  -4.88600 35.77200 36.85500
N  -9.14300 42.33400 41.14400
C  -8.11800 42.62800 40.24500
C  -8.09400 44.00500 39.95700
C  -9.10700 44.61100 40.67700
C  -9.83500 43.52300 41.29800
C  -7.10100 44.80200 39.03200
C  -9.33100 46.10500 40.68600
O  -8.59600 46.89800 40.17700
C  -10.57800 46.61100 41.56500
N  -11.40100 41.07400 42.51100
C  -11.80600 42.38000 42.46200
C  -13.34000 42.40200 42.85400
C  -13.69600 41.00800 43.39700
C  -12.37200 40.22400 43.05500
C  -13.75900 43.56700 43.70500
C  -14.84900 40.42700 42.58800
C  -15.68100 39.35600 43.24200
N  -10.03000 38.45700 42.37900
C  -11.16900 37.95500 42.96200
C  -11.08800 36.54300 43.06100
C  -9.91600 36.15600 42.36900
C  -9.38800 37.38300 41.87700
C  -12.09900 35.62000 43.61300
C  -9.03800 35.08300 41.88300
O  -9.11400 33.86600 42.09900
C  -7.91300 35.74600 40.99900
C  -6.65000 35.63000 41.67800
O  -6.37600 36.35300 42.66500
O  -5.96400 34.63500 41.08200
C  -4.53900 34.42700 41.45900
H  -6.43900 42.13800 39.15400
H  -11.53000 44.45500 42.04300
H  -13.01900 38.25900 43.86600
H  -5.46700 37.65700 40.06900
H  -6.10300 39.87100 38.19000
H  -4.30300 40.89300 39.55000
H  -3.75300 39.31500 38.98800
H  -4.26800 39.50000 40.63400
H  -6.11700 37.88500 37.43800
H  -7.61400 37.12100 38.06300
H  -6.50000 35.14700 38.15800
H  -5.08700 35.86800 38.94600
H  -7.73100 45.35500 38.33500
H  -6.44100 44.13800 38.47400
H  -6.37700 45.43500 39.54500
H  -10.52000 47.69900 41.55500
H  -10.39300 46.15100 42.53600
H  -11.42100 46.18000 41.02500
H  -13.83700 42.54600 41.89400
H  -13.95800 40.94300 44.45300
H  -14.59900 43.25300 44.32400
H  -14.03900 44.35800 43.00900
H  -12.95300 44.01500 44.28500
H  -14.53800 40.03800 41.61800
H  -15.54600 41.19600 42.25500
H  -15.47900 39.25300 44.30800
H  -15.45200 38.41100 42.74900
H  -16.73600 39.48600 43.00100
H  -12.69800 36.09600 44.39000
H  -11.53400 34.78400 44.02600
H  -12.74800 35.36800 42.77400
H  -7.89000 35.41100 39.96200
H  -4.47300 33.34300 41.35500
H  -4.28400 34.79600 42.45200
H  -3.86400 34.76500 40.67200


