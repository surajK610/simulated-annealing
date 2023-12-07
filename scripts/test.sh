module load cmake
module load googletest
module load clang

cd build  
cmake ..  
make     
make test