%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 17.10100 -2.07700 27.77000
C  16.45500 -0.19100 30.58100
C  19.23600 -3.94700 29.68100
C  17.96100 -3.74600 24.96300
C  15.37700 0.23000 25.82200
N  17.83400 -1.91300 29.91400
C  17.31200 -1.24300 30.86800
C  17.86400 -1.55200 32.35400
C  18.99500 -2.58100 31.95700
C  18.76000 -2.80300 30.38800
C  20.50900 -2.04800 32.15900
C  16.79200 -2.13400 33.37400
C  16.96200 -1.78000 34.89400
H  18.25700 -1.14000 35.29400
N  18.31900 -3.63600 27.38700
C  19.10500 -4.38300 28.35000
C  19.74900 -5.47800 27.63200
C  19.48000 -5.36700 26.24500
C  18.51900 -4.24500 26.14900
C  20.50800 -6.46400 28.44500
C  19.91900 -6.24400 25.04900
O  19.71600 -6.00400 23.89600
C  20.84100 -7.42500 25.41400
N  16.84100 -1.78700 25.63700
C  17.26600 -2.56700 24.69400
C  17.04700 -1.94200 23.28700
C  16.12300 -0.73100 23.56000
C  16.07300 -0.74300 25.10600
C  18.39300 -1.56900 22.63000
C  14.74200 -0.86000 22.84700
C  14.78100 -0.68700 21.32000
N  16.09100 -0.30700 28.09200
C  15.35600 0.45700 27.20300
C  14.73300 1.49900 27.97600
C  15.14100 1.35600 29.33600
C  15.94900 0.18700 29.29900
C  13.91400 2.52600 27.44500
C  15.02300 1.80700 30.68600
O  14.43000 2.83800 31.12300
C  15.98100 0.86900 31.51800
C  16.98700 1.71200 32.12800
O  17.99700 2.11700 31.55300
O  16.51500 2.07400 33.36400
C  17.31500 3.14200 34.00000
H  19.94300 -4.54500 30.25900
H  18.26400 -4.34800 24.10300
H  14.78900 0.88500 25.17500
H  18.31000 -0.60500 32.65900
H  18.90200 -3.48800 32.55400
H  20.47800 -0.95900 32.20400
H  21.12500 -2.36600 31.31800
H  21.03200 -2.40200 33.04800
H  16.90400 -3.21300 33.27100
H  15.74800 -2.07700 33.06400
H  16.83100 -2.63200 35.56100
H  16.22600 -1.02300 35.16500
H  21.52400 -6.11200 28.62600
H  20.42400 -7.48300 28.06700
H  20.00400 -6.52900 29.40900
H  20.33800 -8.26100 25.90100
H  21.56700 -6.93500 26.06200
H  21.42700 -7.84100 24.59500
H  16.51200 -2.73300 22.76200
H  16.63000 0.18800 23.26600
H  19.22100 -1.38900 23.31600
H  18.44300 -0.64500 22.05300
H  18.63700 -2.41000 21.98100
H  14.08400 -0.07400 23.21500
H  14.28900 -1.82100 23.09300
H  15.82700 -0.60000 21.02600
H  14.42000 0.31100 21.07100
H  14.22200 -1.38300 20.69500
H  13.25200 2.01200 26.74800
H  14.42200 3.30300 26.87400
H  13.44500 3.05700 28.27400
H  15.32300 0.48100 32.29600
H  17.17000 4.03200 33.38700
H  18.34600 2.83300 34.16900
H  16.87700 3.33200 34.98000
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


