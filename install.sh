chmod 764 */*/*/*
chmod 764 */*/*
chmod 764 */*
chmod 764 *
cd  idba
./build.sh
cd ..
cd metaSpades/bin/
p=$(pwd)
ln -s "$p"/spades.py "$p"/metaspades.py
cd ../../
