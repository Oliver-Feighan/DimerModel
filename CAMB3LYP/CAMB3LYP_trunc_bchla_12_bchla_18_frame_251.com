%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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


