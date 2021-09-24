%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg -1.96800 16.99100 27.07800
C  -2.18500 14.77300 29.87900
C  -3.22300 19.43300 29.10400
C  -2.21300 18.84100 24.51400
C  -1.82200 14.16700 25.04100
N  -2.58500 16.98200 29.20400
C  -2.56200 16.08300 30.18300
C  -2.90400 16.68200 31.58000
C  -3.74100 17.95500 31.11000
C  -3.12300 18.15900 29.73800
C  -5.23100 17.67600 30.94100
C  -1.62300 17.10100 32.37700
C  -1.84400 17.38800 33.86500
H  -0.80000 17.02400 34.81100
N  -2.15300 18.95200 26.92600
C  -2.67100 19.79100 27.85500
C  -2.58800 21.14700 27.27100
C  -2.15500 21.00600 25.85300
C  -2.18200 19.57400 25.69800
C  -2.85500 22.39700 28.05200
C  -1.78700 22.00400 24.80100
O  -1.66500 21.68100 23.57000
C  -1.65800 23.40500 25.24100
N  -2.23300 16.57800 25.04900
C  -2.24700 17.56000 24.14100
C  -2.16400 17.01500 22.74600
C  -1.84400 15.45800 22.89800
C  -1.88200 15.35900 24.44300
C  -3.41900 17.19600 21.88100
C  -0.59500 14.92100 22.15500
C  0.65800 15.50100 22.77200
N  -1.81800 14.84700 27.38100
C  -1.79500 13.85500 26.39000
C  -1.79500 12.55300 27.08200
C  -1.83500 12.88000 28.43600
C  -1.88500 14.29600 28.58000
C  -1.81900 11.27600 26.45500
C  -1.87200 12.36200 29.79400
O  -1.76000 11.24400 30.29900
C  -2.22300 13.56300 30.77900
C  -1.25100 13.50900 31.92600
O  -0.03700 13.73600 31.74700
O  -1.93000 13.43500 33.08100
C  -1.02500 13.40800 34.26100
H  -3.82700 20.21100 29.57500
H  -2.18600 19.44400 23.60400
H  -1.69200 13.24400 24.47300
H  -3.52200 16.04200 32.21000
H  -3.55300 18.82200 31.74300
H  -5.45200 16.61400 31.04700
H  -5.47900 18.01800 29.93600
H  -5.82000 18.29200 31.62100
H  -1.29900 18.08700 32.04400
H  -0.83800 16.35200 32.26700
H  -2.81200 17.01700 34.20300
H  -1.90600 18.46800 33.99900
H  -3.31200 22.17700 29.01600
H  -3.61900 22.95400 27.51000
H  -1.96100 23.01100 28.16200
H  -2.60300 23.73500 25.67300
H  -1.30300 23.96100 24.37400
H  -0.92000 23.45700 26.04200
H  -1.36000 17.49200 22.18700
H  -2.66000 14.83500 22.53100
H  -3.23500 17.73100 20.94900
H  -4.17200 17.70500 22.48200
H  -3.93500 16.27100 21.62500
H  -0.61200 15.23200 21.11100
H  -0.51600 13.83400 22.17600
H  0.74600 16.46900 22.28000
H  1.51900 14.87000 22.54900
H  0.52700 15.68800 23.83800
H  -1.10100 10.65800 26.99400
H  -1.56400 11.34900 25.39800
H  -2.85900 10.99200 26.61500
H  -3.25100 13.41000 31.10600
H  -1.27000 12.51300 34.83300
H  -1.13500 14.37200 34.75800
H  0.03200 13.23700 34.05400
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


