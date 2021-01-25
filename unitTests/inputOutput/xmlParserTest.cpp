#include <catch2/catch.hpp>
#include "InputSchemaEVAA.h"

void parseFile(void) {
	const auto parsedInputFile = InputFileEVAA("inputFileLinear.xml");
}


TEST_CASE("Parse input file", "[xml]") {
    REQUIRE_NOTHROW( parseFile() );
}