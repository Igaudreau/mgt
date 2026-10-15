sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
rm *.eps
rm *.xq
makemg mpi
mg_sn1a < sod.dat
gnuplot plotShockTube.plt
gnuplot plotShockTube_pres.plt
gnuplot plotShockTube_u.plt
xdg-open shockTube.eps
xdg-open shockTube_pres.eps
xdg-open shockTube_u.eps
