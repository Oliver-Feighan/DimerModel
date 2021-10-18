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
Mg -5.36100 24.76000 27.02500
C  -3.64200 26.63800 29.43600
C  -6.40000 22.65900 29.54700
C  -6.94500 22.99100 24.73300
C  -4.07500 26.91800 24.56900
N  -5.09300 24.76200 29.22400
C  -4.33400 25.60400 30.03700
C  -4.64200 25.43000 31.57600
C  -5.37000 24.10100 31.51700
C  -5.66000 23.80300 29.98800
C  -4.43700 22.92100 32.05400
C  -5.65800 26.61300 31.90200
C  -5.45100 27.38800 33.27700
H  -4.92400 26.43700 34.32100
N  -6.35800 23.02700 27.11800
C  -6.68000 22.31300 28.19300
C  -7.48400 21.16300 27.87400
C  -7.81800 21.31800 26.51400
C  -6.97100 22.47200 26.02900
C  -7.92400 20.11100 28.82500
C  -8.85200 20.60400 25.66800
O  -9.21200 21.01500 24.55900
C  -9.41800 19.25800 26.12900
N  -5.43800 24.92500 24.91100
C  -6.26300 24.05500 24.23800
C  -6.35200 24.50400 22.78100
C  -5.50900 25.87000 22.72900
C  -5.01000 25.97700 24.14600
C  -5.90500 23.33400 21.83500
C  -6.23300 27.10700 22.22100
C  -5.47900 27.85000 21.07400
N  -4.10000 26.41500 26.99800
C  -3.63700 27.19900 25.89500
C  -2.75500 28.21500 26.34700
C  -2.74300 28.01700 27.73800
C  -3.53100 26.89000 28.05100
C  -2.03600 29.28900 25.59600
C  -2.22900 28.50400 28.97800
O  -1.46700 29.40000 29.25600
C  -2.89200 27.68300 30.12000
C  -1.83900 27.24500 31.07200
O  -0.96900 26.39800 30.89000
O  -2.04600 27.91300 32.21100
C  -1.18400 27.58800 33.32800
H  -6.90800 22.00000 30.25400
H  -7.48200 22.35700 24.02300
H  -3.70300 27.64900 23.84800
H  -3.68500 25.58900 32.07200
H  -6.31300 24.16400 32.05900
H  -3.44500 23.27200 32.33900
H  -4.33700 22.17800 31.26300
H  -4.99100 22.51200 32.89900
H  -6.69800 26.29000 31.93800
H  -5.59700 27.30300 31.06000
H  -6.41300 27.65600 33.71400
H  -4.98200 28.36000 33.12200
H  -9.00400 20.22200 28.93100
H  -7.47900 20.25200 29.80900
H  -7.70000 19.10000 28.48600
H  -9.79600 18.88600 25.17700
H  -10.18200 19.34000 26.90300
H  -8.70400 18.49400 26.43600
H  -7.35300 24.74000 22.42000
H  -4.63900 25.71000 22.09200
H  -6.65200 23.29700 21.04200
H  -5.78500 22.37600 22.34100
H  -5.04100 23.62000 21.23400
H  -6.46800 27.83800 22.99400
H  -7.10900 26.70500 21.71000
H  -4.60800 27.27900 20.75200
H  -5.00600 28.76500 21.43000
H  -6.13000 28.01900 20.21600
H  -2.85300 29.97800 25.38000
H  -1.55300 28.80400 24.74800
H  -1.32200 29.82700 26.21800
H  -3.47700 28.43300 30.65400
H  -0.13900 27.68100 33.03100
H  -1.23800 26.56800 33.71000
H  -1.41700 28.31200 34.10900


