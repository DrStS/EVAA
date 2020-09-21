***
# EVAA  
EVAA: Efficient Vehicle dynAmics simulAtor
high performance vehicle dynamics model and simulation for applications 
w.r.t. autonomous driving
***
# Status: C++ Continuous Integration   
The status of continuous integration:
![Configure, Build and UnitTest](https://github.com/DrStS/EVAA/workflows/Configure,%20Build%20and%20UnitTest/badge.svg)  
The status of LaTeX documentation:
![Build LaTeX doc](https://github.com/DrStS/EVAA/workflows/Build%20LaTeX%20doc/badge.svg)  
The status of Doxygen documentation:
![Build Doxygen doc](https://github.com/DrStS/EVAA/workflows/Build%20Doxygen%20doc/badge.svg)  
[https://drsts.github.io/EVAA/](https://drsts.github.io/EVAA/)
***
# Build EVAA  
Use cmake to configure EVAA
Windows 
```console
cmake -G "Visual Studio 16 2019" ..
cmake --build . --config Release
ctest
```
