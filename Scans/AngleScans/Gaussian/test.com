%chk=acetone_td
%nproc=4
#p cam-b3lyp/aug-cc-pvtz td=(nstates=6,root=1)

TDDFT excited states in gaussian

0 1
c      -0.13790764   1.26907390   0.15072390 
c       0.00000000   0.00000000  -0.65623430 
o       0.00000000   0.00000000  -1.85927568 
c       0.13790764  -1.26907390   0.15072390 
h      -0.17840414   2.14171260  -0.49664289 
h       0.69579306   1.37757065   0.84709541 
h      -1.04589755   1.23548202   0.75657857 
h       1.04589755  -1.23548202   0.75657857 
h       0.17840414  -2.14171260  -0.49664289 
h      -0.69579306  -1.37757065   0.84709541
