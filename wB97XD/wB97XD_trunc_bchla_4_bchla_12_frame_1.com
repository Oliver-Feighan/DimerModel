%nproc=24
%mem=175gb
#p wB97XD/Def2SVP td=(nstates=5)

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
Mg 47.67100 15.89000 28.45400
C  45.34400 15.65500 31.04800
C  49.48000 18.07500 30.51500
C  49.68700 16.65900 25.89600
C  45.30600 14.62100 26.30000
N  47.50700 16.71700 30.64000
C  46.40900 16.47700 31.40500
C  46.52000 17.23500 32.78700
C  47.75800 18.22000 32.48000
C  48.27100 17.69900 31.12700
C  47.16600 19.65600 32.25100
C  46.90400 16.16800 33.87700
C  47.20400 16.68600 35.31000
H  48.31600 16.29900 36.23900
N  49.43000 17.01400 28.32000
C  50.06900 17.82700 29.27100
C  51.20000 18.39100 28.66300
C  51.38300 17.85100 27.39700
C  50.14500 17.12800 27.15700
C  52.11500 19.27100 29.45000
C  52.51500 18.07600 26.34100
O  52.57400 17.52300 25.24800
C  53.78800 18.86100 26.66900
N  47.44200 15.76700 26.39600
C  48.48400 16.01700 25.57200
C  48.26000 15.53000 24.12700
C  46.80100 15.00400 24.22400
C  46.51000 15.00800 25.71400
C  48.61100 16.58300 23.02400
C  46.48600 13.65000 23.53100
C  46.96600 12.43800 24.21000
N  45.75600 15.16700 28.59900
C  44.93800 14.69700 27.68600
C  43.69800 14.18000 28.27400
C  43.86300 14.56100 29.64400
C  45.08900 15.18100 29.78600
C  42.77800 13.18000 27.74800
C  43.21100 14.49100 30.94000
O  42.19200 13.93500 31.28200
C  44.14500 15.27600 31.86900
C  43.41900 16.43800 32.55500
O  43.17000 17.49000 32.02600
O  43.04200 16.05500 33.82000
C  42.33900 17.13300 34.55300
H  50.09800 18.64700 31.20900
H  50.41400 16.70200 25.08200
H  44.81200 14.14000 25.45400
H  45.63400 17.83000 33.01000
H  48.52200 18.24600 33.25800
H  47.13800 20.01300 31.22200
H  47.84000 20.42900 32.62100
H  46.16600 19.79300 32.66200
H  47.75800 15.61000 33.49400
H  46.09800 15.43400 33.89800
H  46.24800 16.58700 35.82500
H  47.42100 17.73200 35.09700
H  51.80300 19.48300 30.47300
H  52.34700 20.14300 28.83700
H  53.02500 18.67700 29.53900
H  53.61800 19.92500 26.83400
H  54.61200 18.70300 25.97300
H  54.19200 18.41700 27.57900
H  49.00300 14.73300 24.13500
H  46.20400 15.80500 23.78800
H  49.40000 17.25400 23.36200
H  47.72600 17.18100 22.80500
H  48.98600 16.12600 22.10800
H  46.99700 13.77600 22.57700
H  45.42300 13.59700 23.29500
H  46.25100 12.14600 24.97900
H  47.97200 12.58100 24.60400
H  46.98800 11.57800 23.54000
H  43.37600 12.59600 27.04800
H  41.93300 13.49000 27.13300
H  42.41600 12.61300 28.60500
H  44.52600 14.58100 32.61800
H  41.57800 16.67700 35.18600
H  41.70900 17.79000 33.95300
H  43.02300 17.66800 35.21100


