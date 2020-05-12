#include <benchmark/benchmark.h>
#include <mkl.h>
#include <iostream>
#include "EAAAComputeEngine.h"

using namespace EVAA;

class Setup
{
    Setup()
    {
        // call your setup function
        std::cout << "singleton ctor called only once in the whole program" << std::endl;
        /// create db
        db = MetaDatabase<Constants::floatEVAA>::getDatabase();
    }

public:
    MetaDatabase<Constants::floatEVAA>& db;

    static void PerformSetup(const std::string& loadFilename)
    {
        static Setup setup;
        db.readLoadParameters(loadFilename);
    }
};

static void BM_SomeFunction(benchmark::State& state) {
    Setup::PerformSetup()
        for (auto _ : state) {
            // ...
        }
}