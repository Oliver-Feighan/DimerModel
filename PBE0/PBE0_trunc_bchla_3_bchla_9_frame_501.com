%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 1.29600 7.71300 26.88100
C  1.46900 9.82700 29.67000
C  2.42600 5.06900 28.94700
C  1.47300 5.72500 24.13400
C  0.92500 10.50000 24.87700
N  1.69900 7.49100 29.04900
C  1.51500 8.49600 30.00300
C  1.43800 7.84300 31.43000
C  1.83300 6.38800 31.15200
C  2.05200 6.30500 29.60100
C  3.12000 6.01700 31.93700
C  0.04000 8.06800 31.99300
C  -0.12000 8.13500 33.50300
H  1.05800 7.66300 34.38500
N  1.84800 5.68100 26.53500
C  2.31100 4.77600 27.50900
C  2.63500 3.49600 26.90400
C  2.31300 3.59100 25.49300
C  1.87500 5.04400 25.33500
C  3.03900 2.24400 27.69300
C  2.37500 2.55400 24.38800
O  2.03400 2.84400 23.24900
C  2.64700 1.14500 24.75200
N  1.13700 8.09300 24.78900
C  1.15400 7.08500 23.87300
C  0.77800 7.64200 22.47500
C  0.53200 9.18300 22.70900
C  0.76500 9.29400 24.23800
C  1.89300 7.33600 21.45900
C  -0.89600 9.67800 22.35700
C  -0.86300 10.38200 21.00400
N  1.30900 9.79300 27.10700
C  1.22300 10.82400 26.23200
C  1.50800 12.06900 26.98600
C  1.50300 11.75700 28.31400
C  1.42300 10.34100 28.32600
C  1.55300 13.42200 26.29600
C  1.55800 12.26300 29.69200
O  1.68500 13.41500 30.08900
C  1.44800 11.03500 30.57500
C  2.56000 11.13900 31.53300
O  3.76600 10.89700 31.38900
O  2.04800 11.65700 32.70500
C  2.99200 11.90900 33.77200
H  2.69600 4.26200 29.63100
H  1.45200 4.99400 23.32300
H  0.89300 11.33000 24.16800
H  2.21800 8.40600 31.94400
H  1.06000 5.77700 31.61800
H  4.00400 5.85500 31.32000
H  2.97100 4.99800 32.29300
H  3.38200 6.77100 32.67900
H  -0.56100 7.28300 31.53300
H  -0.37400 8.97400 31.55200
H  -1.02500 7.59600 33.78000
H  -0.31800 9.19800 33.64500
H  3.22100 2.42800 28.75200
H  3.95300 1.82100 27.27400
H  2.28900 1.47200 27.52500
H  1.81900 0.77800 25.35900
H  3.66100 1.02700 25.13400
H  2.61600 0.48900 23.88200
H  -0.10600 7.10300 22.13500
H  1.34500 9.68200 22.18000
H  2.33200 8.19000 20.94400
H  1.31000 6.77100 20.73200
H  2.77000 6.84600 21.88100
H  -1.17200 10.51700 22.99600
H  -1.64500 8.88600 22.34300
H  -1.71200 9.95300 20.47200
H  -0.02300 10.09700 20.37100
H  -0.94100 11.46900 21.03900
H  1.59500 13.38900 25.20700
H  2.38500 14.03800 26.63500
H  0.58900 13.91000 26.44400
H  0.47700 11.06900 31.06800
H  3.81000 12.49700 33.35500
H  3.40600 11.00100 34.21100
H  2.51900 12.68100 34.37800
Mg 35.26200 1.75200 29.96200
C  32.85500 2.58300 32.32000
C  37.64200 1.80100 32.35100
C  37.51600 1.12500 27.56800
C  32.74800 1.88400 27.52800
N  35.16600 2.13300 32.11400
C  34.11900 2.35200 32.87300
C  34.63200 2.45300 34.38300
C  36.03300 2.95200 34.09900
C  36.35200 2.21400 32.76000
C  36.05200 4.49500 34.04900
C  34.37500 1.08800 35.15500
C  34.37600 1.20900 36.68200
H  34.84400 2.50900 37.31300
N  37.24000 1.32700 29.98300
C  38.04600 1.32500 31.11900
C  39.36800 0.74800 30.73100
C  39.37200 0.59900 29.27900
C  38.02100 1.10700 28.92100
C  40.46600 0.47300 31.73100
C  40.47100 0.00700 28.28100
O  40.29400 -0.12300 27.06500
C  41.78300 -0.42500 28.90200
N  35.17700 1.65000 27.89700
C  36.22800 1.39600 27.12200
C  35.89000 1.47100 25.60100
C  34.35600 1.46600 25.66700
C  34.03700 1.61900 27.12100
C  36.49900 2.69400 24.86000
C  33.56700 0.32100 24.88900
C  33.00100 -0.79800 25.70000
N  33.26100 2.34500 29.84600
C  32.36300 2.28700 28.83400
C  31.06300 2.71100 29.21500
C  31.19600 2.84400 30.63300
C  32.52400 2.60100 30.90000
C  29.78600 2.98300 28.46300
C  30.54600 3.13100 31.87700
O  29.40100 3.42500 32.14500
C  31.52900 2.83000 33.05300
C  31.55800 3.97100 33.97300
O  31.92600 5.08300 33.62600
O  30.95500 3.71000 35.23000
C  30.85200 4.85600 36.14700
H  38.39100 1.89200 33.13900
H  38.28200 0.96200 26.80800
H  32.00200 1.96400 26.73500
H  34.02600 3.17400 34.93000
H  36.78800 2.59500 34.79900
H  36.64600 4.74400 34.92900
H  35.05700 4.93600 34.12200
H  36.63900 4.75700 33.16900
H  35.14200 0.35900 34.89200
H  33.36400 0.79400 34.87200
H  34.98900 0.39000 37.05700
H  33.35800 1.00500 37.01300
H  40.05900 0.17200 32.69700
H  41.18400 1.27200 31.91300
H  41.03100 -0.39400 31.38900
H  42.08900 0.44800 29.47800
H  42.58200 -0.61900 28.18700
H  41.63400 -1.33400 29.48500
H  36.25600 0.59300 25.06900
H  34.05000 2.44500 25.30000
H  35.69900 3.06000 24.21600
H  37.40400 2.52000 24.28000
H  36.80800 3.41800 25.61400
H  34.19600 -0.03500 24.07300
H  32.66800 0.75100 24.44800
H  31.91500 -0.71000 25.66100
H  33.14600 -0.73200 26.77900
H  33.24300 -1.79400 25.32900
H  29.91400 2.90200 27.38300
H  29.37400 3.93500 28.79900
H  28.98400 2.25500 28.58700
H  31.23100 1.94400 33.61300
H  30.55000 4.60700 37.16400
H  30.10200 5.48600 35.66800
H  31.81400 5.36300 36.08200


