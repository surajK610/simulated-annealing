module load gcc/10.2
module load googletest
module load clang/18

mkdir build -p
cd build
cmake ..
make