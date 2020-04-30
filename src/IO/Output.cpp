// TODO: Copright check

#include "Output.h"

#include <exception>
#include <iostream>

namespace EVAA {
namespace IO {

void checkFileExists(const std::string& filename) {
    std::ifstream f(filename);

    if (f.good()) {
        std::cout << "Filename exists: " << filename;
        return;
    }

    throw std::runtime_error("File " + filename + " does not exist");
}

}  // namespace IO
}  // namespace EVAA
