#include <cmath>
#include "./app4triqs.hpp"

namespace app4triqs {

  toto &toto::operator+=(toto const &b) {
    this->i += b.i;
    return *this;
  }

  toto toto::operator+(toto const &b) const {
    auto res = *this;
    res += b;
    return res;
  }

  bool toto::operator==(toto const &b) const { return (this->i == b.i); }

  void h5_write(h5::group grp, std::string subgroup_name, toto const &m) {
    grp = subgroup_name.empty() ? grp : grp.create_group(subgroup_name);
    h5_write(grp, "i", m.i);
    h5_write_attribute(grp, "Format", toto::hdf5_format());
  }

  void h5_read(h5::group grp, std::string subgroup_name, toto &m) {
    grp = subgroup_name.empty() ? grp : grp.open_group(subgroup_name);
    int i;
    h5_read(grp, "i", i);
    m = toto(i);
  }

  int chain(int i, int j) {
    int n_digits_j = j > 0 ? (int)log10(j) + 1 : 1;
    return i * int(pow(10, n_digits_j)) + j;
  }

} // namespace app4triqs
