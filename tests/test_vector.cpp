#include <gtest/gtest.h>
#include "../core/cvector.hpp"

TEST(CVectorTest, Addition)
{
    CVector<double> empty;
    CVector<double> empty2;

    CVector<double> add = empty + empty2;
    ASSERT_EQ(add.length(), 0.);
    CVector<double> sub = empty - empty2;
    ASSERT_EQ(sub.length(), 0);
}
TEST(CVectorTest, Multiplication)
{
    CVector<double> base(1, 1, 1);
    ASSERT_EQ(base.length(), sqrt(3));
    CVector<double> twobase = base * 2;
    ASSERT_EQ(twobase.length(), sqrt(12));
}
TEST(CVectorTest, VectorAlter)
{
    CVector<double> base2(1, 0, 1);
    ASSERT_EQ(base2.length(), sqrt(2));
    base2.y = 1;
    base2.z = 0;
    base2.x = 0;
    ASSERT_EQ(base2.length(), sqrt(1));
    ASSERT_EQ(base2.y, 1);
    ASSERT_EQ(base2.x, 0);
    ASSERT_EQ(base2.z, 0);
}
