%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 25.49300 50.21300 27.08700
C  23.37900 51.54000 29.65700
C  27.62800 49.17200 29.69500
C  27.85200 49.59100 24.90700
C  23.46200 51.62200 24.75200
N  25.44300 50.31400 29.39800
C  24.41900 50.80900 30.13400
C  24.66700 50.47600 31.59600
C  26.00400 49.68800 31.68900
C  26.38100 49.68400 30.19500
C  27.10100 50.37400 32.68300
C  23.49400 49.61900 32.23900
C  22.96400 50.09200 33.59300
H  23.92100 50.77900 34.64800
N  27.50700 49.52400 27.28400
C  28.12900 49.15000 28.38200
C  29.45200 48.73800 28.00600
C  29.55900 48.68700 26.59300
C  28.31600 49.30500 26.21200
C  30.44200 48.13300 29.00000
C  30.57100 48.09600 25.62400
O  30.49100 48.23300 24.41700
C  31.66500 47.26500 26.23600
N  25.76300 50.88600 25.17100
C  26.76300 50.29700 24.39500
C  26.51800 50.44800 22.86900
C  25.02700 50.95500 22.90400
C  24.70800 51.17300 24.38300
C  27.54000 51.40800 22.24300
C  23.97400 49.88700 22.34200
C  23.37400 50.12000 20.99600
N  23.75000 51.22900 27.14200
C  23.04900 51.77400 26.08800
C  21.90700 52.57500 26.52100
C  21.98600 52.51200 27.95800
C  23.17000 51.80300 28.27000
C  20.85200 53.08800 25.60900
C  21.31100 52.84500 29.20300
O  20.22900 53.44600 29.33300
C  22.22800 52.23700 30.41000
C  22.63900 53.39700 31.18500
O  23.23000 54.34500 30.71600
O  22.29300 53.27500 32.47200
C  22.54800 54.35700 33.47400
H  28.29000 48.76300 30.46100
H  28.36300 49.05600 24.10300
H  22.79700 51.98300 23.96500
H  24.78800 51.42300 32.12200
H  25.86900 48.70800 32.14700
H  26.61300 51.02100 33.41200
H  27.77100 50.92200 32.02000
H  27.67200 49.59700 33.19100
H  23.69200 48.55100 32.31900
H  22.62100 49.66400 31.58800
H  22.46400 49.28700 34.13200
H  22.16000 50.77100 33.30900
H  30.45600 47.04300 29.00900
H  30.32700 48.56400 29.99500
H  31.42700 48.48300 28.69200
H  31.35000 46.60400 27.04400
H  32.50000 47.83500 26.64100
H  31.97600 46.50800 25.51600
H  26.62800 49.45900 22.42500
H  24.97100 51.86600 22.30800
H  28.36600 51.68000 22.90000
H  27.08200 52.35600 21.96300
H  27.94500 50.95800 21.33600
H  23.09600 49.73600 22.97000
H  24.37200 48.87300 22.33000
H  23.66900 49.32700 20.31000
H  23.80100 51.07900 20.70200
H  22.28600 50.13300 21.06000
H  19.87400 53.08000 26.08900
H  20.94000 52.51200 24.68800
H  21.11300 54.12300 25.38400
H  21.65900 51.56400 31.05000
H  22.17700 55.30500 33.08400
H  23.60800 54.48100 33.69700
H  22.16500 54.13200 34.46900
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


