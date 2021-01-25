#include <gtest/gtest.h>
#include "InputSchemaEVAA.h"

void parseFile(void) {
	const auto parsedInputFile = InputFileEVAA("inputFileLinear.xml");
}

TEST(xmlParser, testForNoException) {
    EXPECT_NO_THROW( parseFile() );
}