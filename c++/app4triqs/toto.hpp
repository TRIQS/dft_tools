#pragma once
#include <cpp2py.hpp>
#include <triqs/gfs.hpp>

namespace app4triqs {

  /**
   * A very useful and important class
   *
   * @note A Useful note
   */
  class toto {

    int i = 0.0;

    public:
    toto()  = default;

    /**
     * Construct from integer
     *
     * @param i_ a scalar
     */
    explicit toto(int i_) : i(i_) {}

    ~toto() = default;

    // Copy/Move construction
    toto(toto const &) = default;
    toto(toto &&)      = default;

    /// Copy/Move assignment
    toto &operator=(toto const &) = default;
    toto &operator=(toto &&) = default;

    /// Simple accessor
    int get_i() const { return i; }

    /// Arithmetic operations
    toto operator+(toto const &b) const;
    toto &operator+=(toto const &b);

    /// Comparison
    bool operator==(toto const &b) const;

    /// HDF5
    static std::string hdf5_scheme() { return "Toto"; }

    friend void h5_write(triqs::h5::group grp, std::string subgroup_name, toto const &m);
    friend void h5_read(triqs::h5::group grp, std::string subgroup_name, toto &m);

    /// Serialization
    template <class Archive> void serialize(Archive &ar, const unsigned int version) { ar &i; }
  };

  /**
   * Chain digits of two integers
   *
   * Chain the decimal digits of two integers i and j, and return the result
   *
   * @param :math:`i` The first integer
   * @param :math:`j` The second integer
   * @return An integer containing the digits of both i and j
   *
   * @remark
   */
  int chain(int i, int j);

} // namespace app4triqs
