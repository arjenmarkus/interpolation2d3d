gfortran -c ../src/interpolate_frame.f90
gfortran -c ../src/srfpak.f90
gfortran -c ../src/sphpak.f90
gfortran -c ../src/qshep2d.f90
gfortran -c ../src/qshep3d.f90
gfortran -c ../src/alg623.f90
gfortran -c ../src/alg624.f90
gfortran -c ../src/alg660.f90
gfortran -c ../src/alg661.f90

gfortran -o test_srfpak.exe test_srfpak.f90 interpolate_frame.o srfpak.o alg624.o
gfortran -o test_sphpak.exe test_sphpak.f90 interpolate_frame.o sphpak.o alg623.o
gfortran -o test_qshep2d.exe test_qshep2d.f90 interpolate_frame.o qshep2d.o alg660.o
gfortran -o test_qshep3d.exe test_qshep3d.f90 interpolate_frame.o qshep3d.o alg661.o
