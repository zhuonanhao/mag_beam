# A fast computational software for nonlinnear dynamics of magnetoelastic flexible structures

To run this code, you should have a Linux Ubuntu system

0. Update
sudo apt-get update

1. Install lapack
sudo apt-get install libblas-dev liblapack-dev

2. Install gfortran
sudo apt-get install gfortran

3. Install OpenGL
sudo apt-get install freeglut3-dev

4. Download eigen3
https://gitlab.com/libeigen/eigen/-/releases/3.4.0

5. Move eigen3 to /usr/local/include
sudo mv 'your eigen' /usr/local/include

6. Install mkl

Installation:
1. Use following commands
# Add package repository
sudo apt-get install -y gpg-agent wget
wget -qO - https://repositories.intel.com/graphics/intel-graphics.key | sudo apt-key add -
sudo apt-add-repository 'deb [arch=amd64] https://repositories.intel.com/graphics/ubuntu focal main'

# Install run-time packages
sudo apt-get update
sudo apt-get install intel-opencl-icd intel-level-zero-gpu level-zero intel-media-va-driver-non-free libmfx1

# OPTIONAL: Install developer packages
sudo apt-get install libigc-dev intel-igc-cm libigdfcl-dev libigfxcmrt-dev level-zero-dev
stat -c "%G" /dev/dri/render*
... # for GPU (not required)

cd /tmp
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB

echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update

sudo apt install intel-basekit
sudo apt install intel-hpckit(not required)

Must run this cmd everytime for a new terminal
source /opt/intel/oneapi/setvars.sh


7. Make
g++ -I /usr/local/include/eigen-3.3.7/ main.cpp world.cpp setInput.cpp timeStepper.cpp inertialForce.cpp externalGravityForce.cpp dampingForce.cpp elasticStretchingForce.cpp elasticBendingForce.cpp externalMagneticForce.cpp elasticPlate.cpp -lGL -lglut -lGLU -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -llapack -lgfortran -fopenmp -lpthread -lm -Ofast -o simDER

8. Run 
./simDER option.txt
