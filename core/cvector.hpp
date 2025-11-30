#ifndef CORE_CVECTOR_HPP_
#define CORE_CVECTOR_HPP_

#include <functional> // for function
#include <iostream>   // for operator<<, ostream
#include <sstream>    // for char_traits, basic_stringstream, basic_os..
#include <stdexcept>  // for runtime_error
#include <vector>     // for allocator, vector

/// @brief A simple enum to represent the axis
enum Axis { xaxis, yaxis, zaxis, all, none };

template <typename T> class CVector {

public:
  T x, y, z;
  CVector() {
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
  }

  explicit CVector(const std::vector<T> &vec) {
    if (vec.size() != 3) {
      throw std::runtime_error(
          "Failed to create vector -- passed list was not of len 3!");
    }
    this->x = vec[0];
    this->y = vec[1];
    this->z = vec[2];
  }

  CVector(T x, T y, T z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  CVector(const CVector &v) {
    this->x = v.x;
    this->y = v.y;
    this->z = v.z;
  }

  explicit CVector(const std::function<T()> &generator) {
    // the noise should be independent in each direction
    this->x = generator();
    this->y = generator();
    this->z = generator();
  }

  CVector &operator+=(const CVector &v) {
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
  }

  CVector &operator-=(const CVector &v) {
    this->x -= v.x;
    this->y -= v.y;
    this->z -= v.z;
    return *this;
  }

  CVector operator+(CVector v) {
    CVector res(x + v.x, y + v.y, z + v.z);

    return res;
  };

  CVector operator+(const CVector &v) const {
    CVector res(x + v.x, y + v.y, z + v.z);

    return res;
  };

  CVector operator+(const T &val) const {
    CVector res(x + val, y + val, z + val);
    return res;
  }

  CVector operator-(CVector v) {
    CVector res(x - v.x, y - v.y, z - v.z);

    return res;
  };
  CVector operator-(const CVector &v) const {
    CVector res(x - v.x, y - v.y, z - v.z);

    return res;
  };

  void operator=(CVector v) {
    x = v.x;
    y = v.y;
    z = v.z;
  }

  bool operator==(const CVector &v) {
    if ((x == v.x) && (y == v.y) && (z == v.z))
      return true;
    return false;
  };

  bool operator==(const CVector &v) const {
    if ((x == v.x) && (y == v.y) && (z == v.z))
      return true;
    return false;
  };

  bool operator!=(const CVector &v) {
    if ((x == v.x) && (y == v.y) && (z == v.z))
      return false;
    return true;
  };

  bool operator!=(const CVector &v) const {
    if ((x == v.x) && (y == v.y) && (z == v.z))
      return false;
    return true;
  };

  CVector operator*(const T &val) {
    CVector res(x * val, y * val, z * val);
    return res;
  };

  CVector operator*(const T &val) const {
    const CVector res(x * val, y * val, z * val);
    return res;
  }

  friend CVector operator*(const T &val, const CVector &v) {
    return CVector(val * v.x, val * v.y, val * v.z);
  }

  CVector &operator*=(T v) {
    x *= v;
    y *= v;
    z *= v;
    return *this;
  }

  CVector operator/(T val) {
    if (val == 0) {
      throw std::runtime_error("Failed to divide vector by zero!");
    }
    CVector res(x / val, y / val, z / val);
    return res;
  };

  CVector operator/(T val) const {
    if (val == 0) {
      throw std::runtime_error("Failed to divide vector by zero!");
    }
    CVector res(x / val, y / val, z / val);
    return res;
  };

  T operator[](const int &i) {
    if (i == 0)
      return x;
    else if (i == 1)
      return y;
    else
      return z;
  }

  T operator[](const int &i) const {
    if (i == 0)
      return x;
    else if (i == 1)
      return y;
    else
      return z;
  }

  T length() { return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)); }; // Magnitude

  T length() const {
    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  }; // Magnitude

  void normalize() {
    const T mag = this->length();
    if (mag != 0) {
      x = x / mag;
      y = y / mag;
      z = z / mag;
    }
  };
  void setX(const T &vx) { this->x = vx; }
  void setY(const T &vy) { this->y = vy; }
  void setZ(const T &vz) { this->z = vz; }

  std::vector<T> tolist() { return {this->x, this->y, this->z}; }

  friend std::ostream &operator<<(std::ostream &o, const CVector<T> &obj) {
    o << obj.toString();
    return o;
  }

  std::string toString() {
    std::stringstream ss;
    ss << "[x:" << this->x << ", y:" << this->y << ", z:" << this->z << "]";
    return ss.str();
  }

  const std::string toString() const {
    std::stringstream ss;
    ss << "[x:" << this->x << ", y:" << this->y << ", z:" << this->z << "]";
    return ss.str();
  }

  static CVector<T> fromSpherical(T theta, T phi, T r = 1.0) {
    return CVector<T>(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi),
                      r * cos(theta));
  }
};

#endif // CORE_CVECTOR_HPP_
