%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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


