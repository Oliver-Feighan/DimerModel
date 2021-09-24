%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 53.43000 23.68000 44.03900
C  50.60500 25.72200 43.59200
C  51.43200 20.92900 43.52700
C  56.21400 21.92200 43.71600
C  55.42100 26.61100 44.29000
N  51.34600 23.39300 43.36500
C  50.33600 24.38200 43.24300
C  48.97800 23.81200 43.21200
C  49.26700 22.25100 43.29300
C  50.77800 22.16700 43.44100
C  48.41400 21.54400 44.37400
C  48.22300 24.06100 41.94200
C  49.03400 24.02900 40.63800
H  48.41800 23.46200 39.30000
N  53.81600 21.66600 43.70400
C  52.86000 20.71900 43.65000
C  53.50000 19.42800 43.60900
C  54.89000 19.67000 43.73700
C  55.04900 21.13300 43.75700
C  52.80700 18.10400 43.43000
C  55.97700 18.59900 43.82300
O  55.69800 17.41700 43.77300
C  57.46200 18.80800 43.85800
N  55.50100 24.20200 43.75000
C  56.49500 23.26200 43.68900
C  57.90500 23.90500 43.90100
C  57.59200 25.42600 43.76600
C  56.09200 25.42600 43.96400
C  58.57700 23.60200 45.26200
C  57.80600 26.06200 42.31200
C  58.53500 27.35500 42.24600
N  53.14400 25.68200 44.06800
C  53.99400 26.67600 44.35900
C  53.22600 27.91400 44.76200
C  51.87400 27.51900 44.50800
C  51.89200 26.20900 43.95100
C  53.71800 29.26500 45.26100
C  50.53800 28.03700 44.44100
O  50.10500 29.15300 44.68100
C  49.65300 26.85900 43.86500
C  49.16300 27.30200 42.50200
O  49.80500 27.46200 41.51700
O  47.78800 27.38800 42.55400
C  47.15600 27.52600 41.23800
H  50.76800 20.06300 43.56200
H  57.13800 21.34600 43.78900
H  56.08600 27.44900 44.50900
H  48.42400 24.23500 44.05000
H  49.15700 21.69500 42.36200
H  48.03700 20.69600 43.80200
H  47.64400 22.25100 44.68300
H  49.00900 21.15600 45.20100
H  47.55300 24.92000 41.98700
H  47.44900 23.32600 41.72000
H  49.94100 23.47800 40.88900
H  49.32900 25.07300 40.53100
H  53.31800 17.49200 42.68600
H  51.75400 18.26600 43.20000
H  52.89400 17.62800 44.40700
H  57.89700 19.23300 42.95400
H  57.81600 17.77700 43.85100
H  57.81900 19.19500 44.81300
H  58.56100 23.52900 43.11600
H  58.08000 26.03800 44.52500
H  58.61900 24.52500 45.84000
H  59.48300 23.00800 45.14300
H  57.92700 22.93400 45.82600
H  56.85800 26.15700 41.78300
H  58.38300 25.30400 41.78200
H  59.40900 27.31900 41.59600
H  58.83900 27.73400 43.22200
H  57.92000 28.14400 41.81400
H  53.34500 30.00700 44.55500
H  54.79200 29.18400 45.43300
H  53.17000 29.44900 46.18500
H  48.82700 26.60000 44.52700
H  47.72400 27.04800 40.44000
H  47.12200 28.60500 41.08700
H  46.11800 27.22500 41.37700
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


