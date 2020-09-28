***
# EVAA  
EVAA: Efficient Vehicle dynAmics simulAtor
high performance vehicle dynamics model and simulation for applications 
w.r.t. autonomous driving
***
# Status: C++ Continuous Integration   
The status of continuous integration:
![Configure, Build, UnitTest and Build Doc](https://github.com/DrStS/EVAA/workflows/Configure,%20Build,%20UnitTest%20and%20Build%20Doc/badge.svg)  
***
# Build EVAA  
Install cmake and conan.    
Use cmake to configure EVAA.   
Windows  
```console
cd bin
cmake -G "Visual Studio 16 2019" ..
cmake --build . --config Release
ctest
```
Mac  
```console
cd bin
cmake -G "Xcode" ..
cmake --build . --config Release
ctest
```
# Developer Doc
[https://drsts.github.io/EVAA](https://drsts.github.io/EVAA)
