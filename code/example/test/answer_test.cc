#include <gtest/gtest.h>

#include "allscale_ipic3d/example/answer.h"

using namespace allscale_ipic3d::example;

TEST(AnswerTest, Basic) {
	ASSERT_EQ(42, answer());
}
