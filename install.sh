wd=$(pwd)
chmod 764 */*/*/*
chmod 764 */*/*
chmod 764 */*
chmod 764 *

cd  idba
./build.sh

cd $wd
cd metaspades/bin/
p=$(pwd)
ln -s "$p"/spades.py "$p"/metaspades.py

cd "$wd"/trnascanSE
tar -xvf infernal-1.1.3.tar.gz
cd infernal-1.1.3
./configure
make
make install

cd "$wd"/trnascanSE
tar -xvf automake-1.13.4.tar.gz
cd automake-1.13.4
./configure --prefix=$(pwd)/exec
make
make install

cd "$wd"/trnascanSE
tar -xvf trnascan-se-2.0.6.tar.gz
cd tRNAscan-SE-2.0
./configure --prefix=$(pwd)/exec
cd "$wd"/trnascanSE
eval "sed -i -e 's%missing\ automake\ 1\.13%missing $(pwd)\/automake\-1\.13\.4\/bin\/%g'" ./tRNAscan-SE-2.0/Makefile
cd tRNAscan-SE-2.0
make 
make install

cd "$wd"/trnascanSE
ln -s "$wd"/trnascanSE/infernal-1.1.3/src/* "$wd"/trnascanSE/tRNAscan-SE-2.0/exec/bin/

cd "$wd"/mitfi
tar -xvf infernal-1.0.2.tar.gz
cd infernal-1.0.2
./configure  --prefix=$(pwd)/exec
make
make install
