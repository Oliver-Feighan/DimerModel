%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 8.82800 2.97900 28.39100
C  10.36500 1.60400 31.10700
C  7.59700 5.54500 30.49200
C  7.11400 4.23300 25.80200
C  10.13400 0.47400 26.43700
N  8.69300 3.35100 30.57800
C  9.52500 2.68500 31.48700
C  9.40600 3.27200 32.92200
C  8.46000 4.46400 32.61800
C  8.16100 4.43600 31.13000
C  7.12300 4.69600 33.44300
C  10.76000 3.72700 33.62500
C  10.50900 4.32800 35.03400
H  10.96800 3.41400 36.17100
N  7.54400 4.72800 28.17300
C  7.38800 5.73600 29.10000
C  6.47100 6.69300 28.49200
C  6.33200 6.29300 27.12400
C  7.02200 5.05300 26.95500
C  5.95900 7.85000 29.26000
C  5.41800 7.05100 26.09900
O  5.12000 6.60100 24.97400
C  4.78500 8.36300 26.48300
N  8.63500 2.44800 26.51900
C  7.83100 3.12600 25.56000
C  7.67200 2.29800 24.30500
C  8.94400 1.45600 24.38400
C  9.32300 1.47600 25.85800
C  6.36700 1.43300 24.30400
C  10.11600 2.06600 23.58300
C  10.22200 1.75100 22.11500
N  9.99100 1.29200 28.61900
C  10.57300 0.45400 27.74600
C  11.35400 -0.51100 28.46200
C  11.33200 -0.10700 29.78000
C  10.51600 1.03400 29.82700
C  12.08500 -1.61500 27.78100
C  11.82500 -0.42000 31.05300
O  12.58200 -1.29100 31.43600
C  11.20500 0.70900 32.00300
C  10.44000 0.02800 33.08200
O  9.42100 -0.60800 32.96400
O  10.93800 0.32100 34.28700
C  10.21200 -0.27100 35.47600
H  7.38800 6.37900 31.16500
H  6.56600 4.63400 24.94700
H  10.68600 -0.27500 25.86600
H  8.89500 2.54300 33.55100
H  9.00500 5.40000 32.74400
H  6.85500 5.73300 33.64700
H  7.10100 4.02900 34.30500
H  6.33500 4.37900 32.76000
H  11.20100 4.41000 32.90000
H  11.41200 2.85700 33.70200
H  9.47300 4.65600 35.11900
H  11.13500 5.21700 35.10500
H  4.89200 7.69700 29.42300
H  6.27900 8.73500 28.71000
H  6.46700 7.93200 30.22100
H  4.23800 8.77300 25.63400
H  5.63600 9.03200 26.61800
H  4.16200 8.37000 27.37700
H  7.64200 2.93000 23.41700
H  8.80000 0.42300 24.06800
H  6.45900 0.58300 24.97900
H  6.13800 1.09500 23.29300
H  5.57600 2.10700 24.63100
H  11.08800 1.85700 24.02800
H  10.10000 3.15100 23.68700
H  10.17300 2.70400 21.58900
H  9.27100 1.26300 21.90100
H  11.07000 1.12100 21.84600
H  13.12700 -1.33900 27.62400
H  11.63900 -2.06200 26.89200
H  12.10800 -2.42700 28.50800
H  12.02700 1.28100 32.43400
H  9.15600 -0.02800 35.35600
H  10.52300 0.16000 36.42800
H  10.37400 -1.34900 35.46000
Mg 40.74100 8.40300 29.69900
C  42.66000 10.02500 32.01600
C  38.73200 7.21500 32.07700
C  39.17900 6.75200 27.36300
C  43.28300 9.22800 27.34100
N  40.79300 8.45600 31.86800
C  41.50500 9.39200 32.54400
C  40.94800 9.63500 33.92000
C  39.85300 8.49700 34.00200
C  39.82700 7.94800 32.59200
C  40.10000 7.44100 35.10600
C  40.44100 11.18100 34.18500
C  40.87000 11.88200 35.48300
H  41.23900 10.96600 36.69200
N  39.13800 7.34400 29.69200
C  38.39900 6.91700 30.74900
C  37.34700 6.00300 30.37600
C  37.33700 6.02300 29.01400
C  38.65500 6.62800 28.64700
C  36.61000 5.29800 31.50000
C  36.33400 5.36600 28.04500
O  36.44200 5.38700 26.84100
C  35.20700 4.65100 28.73000
N  41.30000 7.91900 27.67200
C  40.36700 7.31400 26.91500
C  40.80200 7.26400 25.44700
C  41.95100 8.32700 25.38700
C  42.23200 8.49800 26.89400
C  41.21600 5.86100 24.82900
C  41.58000 9.72400 24.76100
C  40.48400 10.59900 25.31000
N  42.57200 9.48100 29.57000
C  43.55300 9.66900 28.59100
C  44.67900 10.46300 29.09900
C  44.36400 10.63700 30.45500
C  43.08700 10.02200 30.64800
C  45.85600 11.01300 28.33500
C  44.78800 11.17000 31.66800
O  45.80500 11.75000 32.01500
C  43.61600 11.01200 32.67100
C  44.29700 10.49000 33.87800
O  44.63700 9.31800 34.05700
O  44.61300 11.49900 34.74100
C  45.54500 11.24000 35.84300
H  37.99000 7.00700 32.85100
H  38.68100 6.27000 26.52000
H  43.95700 9.44200 26.50800
H  41.68000 9.37000 34.68300
H  38.87700 8.95900 34.14900
H  40.48500 6.53200 34.64400
H  39.14600 7.20600 35.57800
H  40.83000 7.87400 35.79000
H  39.35300 11.12200 34.15300
H  40.80700 11.84200 33.40000
H  40.04300 12.46600 35.88600
H  41.65400 12.60800 35.26800
H  37.22900 4.66200 32.13300
H  35.89000 4.51000 31.28000
H  36.12800 6.14400 31.99000
H  34.67500 5.32500 29.40100
H  35.72300 3.80100 29.17600
H  34.37800 4.27900 28.12900
H  39.88500 7.56300 24.93800
H  42.83200 7.86300 24.94300
H  41.11200 5.08600 25.58800
H  42.28100 5.83900 24.60100
H  40.70000 5.56000 23.91600
H  41.36300 9.40900 23.74000
H  42.43100 10.40400 24.73900
H  39.92000 10.04700 26.06300
H  39.77700 10.86400 24.52400
H  40.79900 11.50000 25.83600
H  46.78200 10.48700 28.56300
H  45.97200 12.07500 28.55200
H  45.60900 10.77100 27.30100
H  43.07400 11.94700 32.81400
H  46.05900 12.14500 36.16900
H  46.18500 10.38200 35.64300
H  45.04400 10.86200 36.73400


