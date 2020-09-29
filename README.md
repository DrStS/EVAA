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

# Naming conventions C++: CamelCase
A snippet taken from this [source](https://lefticus.gitbooks.io/cpp-best-practices/content/03-Style.html) mainly.   
* Types start with upper case: <code>MyClass</code>  
* Files start with upper case: <code>MyFile</code>  
* Functions and variables start with lower case: <code>myMethod</code>   
* Constants are all upper case: <code>const double PI=3.14159265358979323;</code>  
* Macro names use upper case with underscores: <code>INT_MAX</code>  
* Template parameter names use camel case: <code>InputIterator</code>  
* All other names use snake case: <code>unordered_map</code>   
Distinguish Private Object Data   
* Name private data with a m_ prefix to distinguish it from public data. m_ stands for "member" data.  
Distinguish Function Parameters   
* Name function parameters with an t_ prefix. t_ can be thought of as "the", but the meaning is arbitrary.  
## Example
Header file MyClass.h    
```cpp
class MyClass
{
public:
  MyClass(MyOtherClass t_myOtherClass)
    : m_myOtherClass(t_myOtherClass)
  {
  }

private:
  MyOtherClass m_myOtherClass;
};
```
Implementation file MyClass.cpp    
```cpp
class MyClass
{
public:
  MyClass(MyOtherClass t_myOtherClass)
    : m_myOtherClass(t_myOtherClass)
  {
  }

private:
  MyOtherClass m_myOtherClass;
};
```
  
