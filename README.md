# rbnics-examples

💾 A personal repository to store reduced basis codes on RBniCS. 


## 🔨 Installation

1. Install [FEniCS](https://fenics.readthedocs.io/en/latest/installation.html#debian-ubuntu-packages) on Ubuntu via Ubuntu Personal Package Archives (PPA)

```bash
sudo
```

```bash
sudo apt-get install --no-install-recommends software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install --no-install-recommends fenics
```

2. Install [PETSc](https://www.mcs.anl.gov/petsc/) and [SLEPc](https://slepc.upv.es/)

* Prerequisites

```bash
sudo apt-get install valgrind
sudo apt-get install gfortran
sudo apt-get install python3-distutils
```

* Download **PETSc** and **SLEPc** (should be in /home/USER/Downloads)
* Install **petsc-3.14.1** locally on a user-defined folder

```bash
mkdir /home/USER/Packages
cd Packages
tar -xzf /home/USER/Downloads/petsc-3.14.1.tar.gz
cd petsc-3.14.1
./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --download-fblaslapack
make all check
```

* Install **slepc-3.14.0** locally on a user-defined folder

```bash
cd /home/USER/Packages
tar -xzf /home/USER/Downloads/slepc-3.14.0.tar.gz
cd slepc-3.14.0
export PETSC_DIR=/home/USER/Packages/petsc-3.14.1
export PETSC_ARCH=arch-linux2-gnu-c-debug
./configure
make all check
```

3. Install **mpi** and **pip3**

```bash
sudo apt install mpi
sudo apt install python3-pip
```

4. Install **petsc4py** and **slepc4py**

```bash
pip3 install petsc4py
pip3 install slepc4py
```

5. Install [RBniCS](https://www.rbnicsproject.org/)

* Go back to **Packages** directory

```bash
git clone https://github.com/RBniCS/RBniCS.git
cd RBniCS
sudo python3 setup.py install
```



<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE2MTIyMzg2NzgsMTE5MzE1MTAxNyw3NT
U1NTEsLTE4ODI3MTgyMzMsLTU5NDkwMDE4N119
-->