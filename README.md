# C++ library for schunk powerball lwa4p (linux)
### contributor: Amir Memar , Sri Sadan 

The code is developed by [HILS lab](https://github.com/ubhils).


# Required libraries
### Hardware
- ntcan (CAN-USB) : https://esd.eu/en 
(Download the latest driver from https://esd.eu/en/software-downloads/27071 - ntcanSck (or) install it from the ESD driver CD). 
  - Socketcan : http://www.pengutronix.de/software/libsocketcan/download/
  - iproute2  : https://github.com/shemminger/iproute2
    - iproute2 might require flex and bison (sudo apt-get install flex; sudo apt-get install bison). "ip link show : should show the connected can device"
    if you face issues with the installation (aligned_u64), try this https://github.com/cloudflare/bpftools/commit/b1bbc6b2a35be84d38355f5e9da1382554ba7f70

- Phidget21 : https://www.phidgets.com/downloads/phidget21/libraries/linux/libphidget/

- if you mistakelnly created the can0 file and want to remove, then open sys/class/net - check if there is any can0 folder if system is not attched to the robot. then you must have to delete it using the following command. 
```
sudo ip link delete dev can0 type vcan
```

- To check the hardware is running correctly, attach robot to the system then run the following command in the bin folder. 


```
sudo ./homing    # if there is no error then it is working perfectly. 
```
### Software
- Boost library : sudo apt-get install libboost-all-dev
- BLAS and LAPACK : sudo apt-get install libblas-dev liblapack-dev
- OpenBLAS : http://www.openblas.net/
- Matplotlib-cpp : https://github.com/lava/matplotlib-cpp
- TooN library : Edward Rosten - https://github.com/edrosten/TooN
- Armadillo : http://arma.sourceforge.net/download.html
- CSV stream : https://github.com/vincentlaucsb/csv-parser#integration
  A useful cpp file for reading the csv files

# Instructions to install the package
```
- git clone https://github.com/srisadhan/powerball.git
- cd powerball
- mkdir build
- cmake CMakeLists.txt -Bbuild
- make Makefile
- make
```

If there is any error of not found sgemm,dgemm, cheev-, then follow this procedure and run previos step again. 

- https://cmake.org/pipermail/cmake/2012-April/049818.html




# Instructions for the Weiss force sensor
#### For using force sensor along with ethernet in linux, modify the /etc/network/interfaces to
```
auto lo
iface lo inet loopback

# Obtain DHCP address from server  
auto en01
iface en01 inet dhcp

# Connect to 192.168.1.10 network over the pcie-ethernet card for force sensor
auto enp2s0
iface enp2s0 inet static
    address 192.168.1.10
    netmask 255.255.255.0
```
replace the en01, enp2s0 to the ethernet adapters recognized on your computer using 'ifconfig'

For more information, please check: https://stackoverflow.com/questions/42922949/two-wired-connection-at-the-same-time


# How to include dynamixel libraries in your own projects
 - Install dynamixel sdk
 - place the "DynamixelSDK/c++/include" folder in your current project or include the directory "/usr/local/include/dynamixel_sdk" in your project
 - Run the following command 
```
g++ read_write.cpp -o read_write -ldxl_x64_cpp -lrt -I include/dynamixel_sdk

(or)

g++ read_write.cpp -o read_write -ldxl_x64_cpp -lrt -I <dynamixel include files>
```

The linking libraries are -lrt (for multi-threading) and -dxl_x64_cpp (dynamixel libraries)




