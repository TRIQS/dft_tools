#include <triqs/test_tools/gfs.hpp>
#include <app4triqs/app4triqs.hpp>

using namespace app4triqs;

TEST(Toto, Add) { // NOLINT

  toto a(0);
  toto b(2);

  auto c = a + b;
  EXPECT_EQ(c, b); // NOLINT
}

TEST(Toto, H5) { // NOLINT

  toto a(0);
  { // Local scope for file
    triqs::h5::file f("f.h5", H5F_ACC_TRUNC);
    h5_write(f, "a", a);
  }

  toto a2;
  {
    triqs::h5::file f("f.h5", H5F_ACC_RDWR);
    h5_read(f, "a", a2);
  }

  EXPECT_EQ(a, a2); // NOLINT
}
