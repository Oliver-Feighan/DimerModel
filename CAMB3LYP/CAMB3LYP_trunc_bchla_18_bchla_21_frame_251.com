%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 34.81300 49.43200 25.14500
C  35.01600 47.65700 28.22400
C  33.32400 52.02800 26.83800
C  34.65200 50.90000 22.18900
C  35.90600 46.46600 23.43400
N  34.14900 49.79800 27.34500
C  34.46100 48.93800 28.38200
C  33.98800 49.65100 29.69500
C  33.38500 51.08500 29.24000
C  33.61700 50.95700 27.73100
C  31.94200 51.39000 29.69300
C  35.01900 49.70600 30.89700
C  34.55900 49.07300 32.17800
H  33.11300 48.96700 32.58500
N  34.17200 51.35300 24.58700
C  33.58900 52.24400 25.43200
C  33.45200 53.46600 24.68200
C  33.72500 53.14400 23.33500
C  34.19500 51.72200 23.26900
C  32.90900 54.73200 25.36600
C  33.46900 53.99900 22.13100
O  33.64100 53.56300 20.96800
C  32.94500 55.37500 22.31400
N  35.18300 48.68800 23.08400
C  35.11600 49.51200 22.09700
C  35.57200 48.95700 20.64400
C  35.95000 47.52900 21.06600
C  35.68700 47.53300 22.55200
C  34.36900 49.06600 19.63700
C  37.40800 47.05400 20.73300
C  38.54700 47.67100 21.57900
N  35.21400 47.49800 25.67400
C  35.61800 46.42200 24.84900
C  35.92400 45.27900 25.72000
C  35.71000 45.72300 27.04700
C  35.31700 47.10200 26.97200
C  36.42700 43.88500 25.23900
C  35.78900 45.30100 28.34400
O  36.05400 44.22400 28.83800
C  35.40700 46.52600 29.22100
C  34.34200 46.13400 30.11500
O  33.13700 46.19100 29.84500
O  34.87800 45.46700 31.22300
C  33.86400 44.85200 32.11500
H  32.63300 52.79500 27.19500
H  34.43900 51.33900 21.21200
H  36.24500 45.54000 22.96500
H  33.26200 48.90600 30.01800
H  33.96500 51.87500 29.71500
H  31.90300 52.32700 30.25000
H  31.43700 50.62000 30.27600
H  31.39300 51.49300 28.75700
H  35.30500 50.73900 31.09600
H  35.87800 49.12400 30.56300
H  34.83800 49.76200 32.97600
H  35.03100 48.09400 32.25300
H  32.96100 54.66600 26.45300
H  31.87600 54.89800 25.06300
H  33.55400 55.53200 25.00400
H  33.66100 55.90400 22.94300
H  31.99400 55.37100 22.84800
H  32.75700 55.84600 21.35000
H  36.44300 49.49500 20.26900
H  35.19100 46.83200 20.71000
H  33.44700 48.97400 20.21100
H  34.41800 48.21400 18.95800
H  34.29100 49.96700 19.02800
H  37.57700 47.10800 19.65800
H  37.43700 45.97400 20.87400
H  38.40300 48.69400 21.92500
H  39.48100 47.54600 21.03000
H  38.64300 47.01100 22.44100
H  37.30600 43.53600 25.78100
H  36.62200 43.84300 24.16800
H  35.74400 43.11800 25.60500
H  36.33300 46.69700 29.76900
H  33.38300 44.00700 31.62300
H  33.04700 45.51600 32.39500
H  34.36200 44.67800 33.06900
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


