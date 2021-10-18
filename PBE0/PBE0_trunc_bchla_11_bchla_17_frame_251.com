%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg 29.10200 58.37000 41.70600
C  26.21100 56.85000 40.52500
C  31.00300 55.92400 39.98500
C  31.78200 60.36200 41.82300
C  27.09800 61.07100 42.76300
N  28.62400 56.74200 40.07100
C  27.36300 56.19900 40.01300
C  27.34200 54.94600 39.19500
C  28.87800 54.50900 39.24100
C  29.59900 55.79900 39.77000
C  29.21600 53.21300 40.13700
C  26.66600 55.21800 37.75400
C  26.06800 54.05700 37.13700
H  26.60500 53.76400 35.72000
N  31.13700 58.12700 41.18100
C  31.72600 57.03700 40.48100
C  33.15300 57.29600 40.39000
C  33.42500 58.57400 40.93600
C  32.09000 59.05300 41.42900
C  34.07200 56.38900 39.66000
C  34.85200 59.26600 40.92400
O  35.77500 58.62800 40.42000
C  35.15900 60.50100 41.76400
N  29.32700 60.43000 42.12600
C  30.57900 60.94600 42.21600
C  30.58100 62.56100 42.56700
C  29.07200 62.72600 43.08400
C  28.42300 61.36300 42.68500
C  31.60200 63.02200 43.67500
C  28.34100 63.90500 42.43900
C  27.75100 65.04300 43.30600
N  27.08200 58.74100 41.92000
C  26.43600 59.84700 42.39100
C  25.02600 59.57700 42.45300
C  24.91200 58.38300 41.64600
C  26.21500 57.96500 41.32600
C  23.95900 60.50800 42.95200
C  23.95200 57.50800 41.00600
O  22.77000 57.42900 41.15100
C  24.75000 56.37400 40.25200
C  24.55300 54.95400 40.82700
O  24.66200 54.66400 41.99900
O  24.13200 54.08300 39.80100
C  23.73600 52.76500 40.19600
H  31.53000 55.10300 39.49600
H  32.55200 61.12000 41.66600
H  26.51800 61.92400 43.12000
H  26.89400 54.10800 39.73000
H  29.20300 54.31600 38.21800
H  29.65600 52.49600 39.44400
H  28.33800 52.79800 40.63100
H  29.96300 53.53300 40.86400
H  27.50600 55.54800 37.14300
H  26.02700 56.07100 37.98400
H  24.99800 54.25700 37.09500
H  26.28800 53.16600 37.72500
H  33.60300 55.42900 39.44500
H  34.92200 56.13300 40.29300
H  34.37800 56.81900 38.70600
H  34.62000 61.39400 41.44600
H  36.23600 60.66600 41.71900
H  34.91600 60.25600 42.79800
H  30.70200 63.02600 41.58900
H  28.95600 62.73400 44.16800
H  32.36800 63.65500 43.22800
H  32.08100 62.14600 44.11400
H  31.02600 63.57600 44.41700
H  27.57000 63.56900 41.74600
H  28.92000 64.44500 41.69000
H  26.66500 65.11700 43.36900
H  28.07200 65.96200 42.81500
H  28.11700 64.98900 44.33200
H  23.07500 59.90100 43.14900
H  23.49000 61.12600 42.18700
H  24.26100 61.15900 43.77200
H  24.32100 56.40300 39.25000
H  23.95500 52.21700 39.27900
H  22.65500 52.67800 40.30400
H  24.24700 52.46200 41.10900


