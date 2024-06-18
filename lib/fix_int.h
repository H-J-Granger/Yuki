#ifndef YUKI_LIB_FIX_INT_H_
#define YUKI_LIB_FIX_INT_H_

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstring>
#include <string>
#include <type_traits>

#include "types.h"

namespace yuki {

/**
 * @brief Handles integer value from -2^len to 2^len.
 *
 * @tparam len Length (in bit), should be a multiple of 64.
 */
template <size_t len>
class fix_int {
  static_assert(len > 0, "The length of fix_int should be greater than 0.");
  static_assert(len % 64 == 0,
                "The length of fix_int should be a multiple of 64.");
  static_assert(len < (1ll << 31),
                "The length of fix_int should be less than 2^31.");

 public:
  static constexpr size_t data_len = len / 64;

  uint64 data[data_len];
  bool sign;  // true if the value is negative. 0 is always non-negative.

 private:
  /**
   * @brief this->data += rhs.data.
   * this->data is supposed to be large enough to store the answer.
   */
  template <size_t len1>
  void add_data(fix_int<len1> const& rhs);

  /**
   * @brief this->data -= rhs.data.
   * Assuming this->compare(rhs) >= 0.
   */
  template <size_t len1>
  void subtract_data(fix_int<len1> const& rhs);

  /**
   * @brief this->data = lhs.data * rhs.data. Using long multiplication.
   * this->data is supposed to be large enough to store the answer.
   */
  template <size_t len1, size_t len2>
  void long_multiplication(fix_int<len1> const& lhs, fix_int<len2> const& rhs);

  /**
   * @brief this->data = lhs.data * rhs.data. Using Karatsuba algorithm.
   * Note that the Karatsuba algorithm is only used once, the multiplication of
   * subroutines are implemented by calling the operator* functions and let it
   * decide. Requires len1 >= 512 && len2 >= 512.
   */
  template <size_t len1, size_t len2>
  void karatsuba(fix_int<len1> const& lhs, fix_int<len2> const& rhs);

  /**
   * @brief this->data = lhs.data * rhs.data. Using Toom-Cook (Toom-3) 
   * algorithm. Note that the Toom-Cook algorithm is only used once, the 
   * multiplication of subroutines are implemented by calling the operator* 
   * functions and let it decide. Requires len1 >= 1024 && len2 >= 1024.
   */
  template <size_t len1, size_t len2>
  void toom_cook(fix_int<len1> const& lhs, fix_int<len2> const& rhs);

  /**
   * @brief this->data = lhs.data * rhs.data. Using FFT over ring F_p, where
   * p = 2524775926340780033 (number_theoretic_transform_t::mod) with a 
   * primitive root of 3. Note that the coefficients should be in [0, 2e18] 
   * and the degree of polynomials should be no greater than 2^24.
   */
  template <size_t len1, size_t len2>
  void number_theoretic_transform(fix_int<len1> const& lhs,
                                  fix_int<len2> const& rhs);

  /**
   * @brief Helper functions and variables for number_theoretic_transform.
   */

  class number_theoretic_transform_t {
   public:
    static constexpr uint64 mod = 2524775926340780033llu;
    /**
     * @brief A primitive root of mod.
     */
    static constexpr uint64 proot = 3;

    /**
     * @brief Montgomery reduction helper class.
     */
    class montgomery_t {
     public:
      static constexpr uint64 mod = number_theoretic_transform_t::mod;
      /**
       * @brief The multiplicative inverse of mod modular 2^64.
       */
      static constexpr uint64 inv_mod = 15348040669855744001llu;
      /**
       * @brief Montgomery reduce.
       */
      static uint64 reduce(uint128 const& x);
      /**
       * @brief Multiply two numbers in Montgomery format.
       */
      static uint64 multiply(uint64 x, uint64 y);
      /**
       * @brief Convert number into Montgomery format.
       */
      static uint64 to_montgomery(uint64 x);
      /**
       * @brief Get number from Montgomery format.
       */
      static uint64 from_montgomery(uint64 x);
    };

    /**
     * @brief the bit inverse of all 15-bit numbers.
     */
    static std::vector<uint64> calc_bit_inv();
    static std::vector<uint64> bit_inv;

    static uint64 bit_inverse(uint64 x);

    /**
     * @brief Modular power.
     */
    static uint64 modular_pow(uint64 x, uint64 y);
    /**
     * @brief Modular inverse.
     */
    static uint64 modular_inv(uint64 x);
    /**
     * @brief Transform into/from NTT format. a.size() must be a power of 2.
     * 
     * @param on false if into NTT format.
     */
    static void transform(std::vector<uint64>& a, bool on);
  };

  /**
   * @brief this->data = lhs.data / rhs.data. Using schoolbook division.
   * this->data is supposed to be large enough to store the answer.
   */
  template <size_t len1, size_t len2>
  void schoolbook_division(fix_int<len1> lhs, fix_int<len2> rhs);

  /**
   * @brief Computes the minimum iteration counts for Newton-Raphson division 
   * process.
   * Floating-point error will more likely to affect on len = 535744, 4286016,
   * 8572032, etc., since the value is pretty close to an integer.
   */
  static constexpr int newton_raphson_times();

 public:
  class fix_point_t {
   public:
    fix_int<len> value;
    int offset;

    fix_point_t();
  };

  template <size_t len1, size_t len2>
  typename fix_int<len1 + len2>::fix_point_t add_fix_point(
      typename fix_int<len1>::fix_point_t const& lhs,
      typename fix_int<len2>::fix_point_t const& rhs);
  template <size_t len1, size_t len2>
  typename fix_int<len1 + len2>::fix_point_t mul_fix_point(
      typename fix_int<len1>::fix_point_t const& lhs,
      typename fix_int<len2>::fix_point_t const& rhs);
  template <size_t len1>
  typename fix_int<len1>::fix_point_t trunc_fix_point(
      typename fix_int<len1>::fix_point_t const&, int off);

  bool get_bit(int index);

 private:
  /**
   * @brief this->data = lhs.data / rhs.data.
   * Using Newton-Raphson division.
   * this->data is supposed to be large enough to store the answer.
   */
  template <size_t len1, size_t len2>
  void newton_raphson_division(fix_int<len1> const& lhs,
                               fix_int<len2> const& rhs);

  /**
   * @brief Only tests data.
   */
  bool is_zero() const;

  static char to_char(int);
  static int from_char(char);

 public:
  fix_int();

  template <std::unsigned_integral T>
  fix_int(T);
  template <std::signed_integral T>
  fix_int(T);

  template <size_t len2>
  fix_int(fix_int<len2> const&);
  template <size_t len2>
  fix_int(fix_int<len2>&&);

  /**
   * @brief Returns the index of the first significant chunk.
   * Returns -1 if data == 0.
   */
  size_t first_significant_chunk() const;

  /**
   * @brief Returns 1 if *this.data > rhs.data, -1 if less, 0 if equal.
   */
  template <size_t len1>
  requires(len >= len1) int compare(fix_int<len1> const& rhs) const;
  // clang-format off
  // this 2 > 1 is for coloring of template stuff.
  template <size_t len1>
  requires(len < len1 && 2 > 1) int compare(fix_int<len1> const& rhs) const;
  // clang-format on

  /**
   * @brief Behaves NOT as built-in integers. The sign is preserved, and the
   * value is shrinked.
   */
  template <size_t len2>
  fix_int<len>& operator=(fix_int<len2> const&);
  template <size_t len2>
  fix_int<len>& operator=(fix_int<len2>&&);

  /**
   * @brief Returns a hex string, without prefix '0x', with lower cases letters.
   */
  std::string to_hex_string();
  /**
   * @brief The letters can be either lower cases or upper cases.
   */
  void from_hex_string(std::string);

  template <size_t len1>
  friend fix_int<len1> operator+(fix_int<len1> const&) noexcept;
  template <size_t len1>
  friend fix_int<len1> operator-(fix_int<len1>);

  template <size_t len1>
  friend fix_int<len1> operator++(fix_int<len1>&);
  template <size_t len1>
  friend fix_int<len1> operator++(fix_int<len1>&, int);
  template <size_t len1>
  friend fix_int<len1> operator--(fix_int<len1>&);
  template <size_t len1>
  friend fix_int<len1> operator--(fix_int<len1>&, int);

  template <size_t len1, size_t len2>
  friend fix_int<len1>& operator+=(fix_int<len1>&, fix_int<len2> const&);
  template <size_t len1, size_t len2>
  friend fix_int<len1>& operator-=(fix_int<len1>&, fix_int<len2> const&);
  template <size_t len1, size_t len2>
  friend fix_int<len1>& operator*=(fix_int<len1>&, fix_int<len2> const&);

  template <size_t len1, size_t len2>
  friend fix_int<len1>& operator/=(fix_int<len1>&, fix_int<len2> const&);
  template <size_t len1, std::unsigned_integral T>
  friend fix_int<len1>& operator/=(fix_int<len1>&, T);
  template <size_t len1, std::signed_integral T>
  friend fix_int<len1>& operator/=(fix_int<len1>&, T);
  template <size_t len1, size_t len2>
  friend fix_int<len1>& operator%=(fix_int<len1>&, fix_int<len2> const&);
  template <size_t len1, std::unsigned_integral T>
  friend fix_int<len1>& operator%=(fix_int<len1>&, T);
  template <size_t len1, std::signed_integral T>
  friend fix_int<len1>& operator%=(fix_int<len1>&, T);

  // operator bool() const;
  // operator uint64() const;

  template <size_t len1, size_t len2>
  friend fix_int<std::max(len1, len2)> operator&(fix_int<len1> const&,
                                                 fix_int<len2> const&);
  template <size_t len1, size_t len2>
  friend fix_int<std::max(len1, len2)> operator|(fix_int<len1> const&,
                                                 fix_int<len2> const&);
  template <size_t len1, size_t len2>
  friend fix_int<std::max(len1, len2)> operator^(fix_int<len1> const&,
                                                 fix_int<len2> const&);
  template <size_t len1>
  friend fix_int<len1> operator~(fix_int<len1>);

  template <size_t len1>
  friend fix_int<len1> operator>>(fix_int<len1> const&, size_t);
  template <size_t len1>
  friend fix_int<len1> operator<<(fix_int<len1> const&, size_t);

  template <size_t len1, size_t len2>
  friend fix_int<len1>& operator&=(fix_int<len1>&, fix_int<len2> const&);
  template <size_t len1, size_t len2>
  friend fix_int<len1>& operator|=(fix_int<len1>&, fix_int<len2> const&);
  template <size_t len1, size_t len2>
  friend fix_int<len1>& operator^=(fix_int<len1>&, fix_int<len2> const&);

  template <size_t len1>
  friend fix_int<len1>& operator>>=(fix_int<len1>&, size_t);
  template <size_t len1>
  friend fix_int<len1>& operator<<=(fix_int<len1>&, size_t);

  template <size_t len1, size_t len2>
  friend fix_int<std::max(len1, len2) + 64> operator+(fix_int<len1> const&,
                                                      fix_int<len2> const&);
  template <size_t len1, size_t len2>
  friend fix_int<std::max(len1, len2) + 64> operator-(fix_int<len1> const&,
                                                      fix_int<len2> const&);
  template <size_t len1, size_t len2>
  friend fix_int<len1 + len2> operator*(fix_int<len1> const&,
                                        fix_int<len2> const&);
  template <size_t len1>
  friend fix_int<len1 + 64> operator*(fix_int<len1>, uint64);

  template <size_t len1, std::unsigned_integral T>
  friend std::pair<fix_int<len1>, fix_int<64>> division(fix_int<len1> const&,
                                                        T);
  template <size_t len1, std::unsigned_integral T>
  friend fix_int<len1> operator/(fix_int<len1> const&, T);
  template <size_t len1, std::unsigned_integral T>
  friend fix_int<64> operator%(fix_int<len1> const&, T);

  /**
   * @brief Divides towards 0.
   */
  template <size_t len1, std::signed_integral T>
  friend std::pair<fix_int<len1>, fix_int<64>> division(fix_int<len1> const&,
                                                        T);
  template <size_t len1, std::signed_integral T>
  friend fix_int<len1> operator/(fix_int<len1> const&, T);
  template <size_t len1, std::signed_integral T>
  friend fix_int<64> operator%(fix_int<len1> const&, T);

  /**
   * @brief Divides towards 0.
   */
  template <size_t len1, size_t len2>
  friend std::pair<fix_int<len1>, fix_int<std::min(len1, len2)>> division(
      fix_int<len1> const&, fix_int<len2> const&);
  template <size_t len1, size_t len2>
  friend fix_int<len1> operator/(fix_int<len1> const&, fix_int<len2> const&);
  template <size_t len1, size_t len2>
  friend fix_int<std::min(len1, len2)> operator%(fix_int<len1> const&,
                                                 fix_int<len2> const&);
};

template <size_t len>
fix_int<len>::fix_int() {
  std::memset(data, 0x00, sizeof(data));
  sign = false;
}

template <size_t len>
template <std::unsigned_integral T>
inline fix_int<len>::fix_int(T x) {
  std::memset(data, 0x00, sizeof(data));
  data[0] = x;
  sign = false;
}
template <size_t len>
template <std::signed_integral T>
inline fix_int<len>::fix_int(T x) {
  std::memset(data, 0x00, sizeof(data));
  if (x < 0) {
    sign = true;
    data[0] = -x;
  } else {
    sign = false;
    data[0] = x;
  }
}

template <size_t len>
template <size_t len2>
inline fix_int<len>::fix_int(fix_int<len2> const& x) {
  sign = x.sign;
  if (len <= len2) {
    std::memcpy(data, x.data, sizeof(data));
  } else {
    std::memcpy(data, x.data, x.data_len * 8);
    std::memset(data + x.data_len, 0x00, 8 * (data_len - x.data_len));
  }
}
template <size_t len>
template <size_t len2>
inline fix_int<len>::fix_int(fix_int<len2>&& x) {
  sign = x.sign;
  if (len <= len2) {
    std::memcpy(data, x.data, sizeof(data));
  } else {
    std::memcpy(data, x.data, x.data_len * 8);
    std::memset(data + x.data_len, 0x00, 8 * (data_len - x.data_len));
  }
}

template <size_t len>
template <size_t len2>
fix_int<len>& fix_int<len>::operator=(fix_int<len2> const& x) {
  sign = x.sign;
  if (len <= len2) {
    std::memcpy(data, x.data, sizeof(data));
  } else {
    std::memcpy(data, x.data, x.data_len * 8);
    std::memset(data + x.data_len, 0x00, 8 * (data_len - x.data_len));
  }
  return *this;
}
template <size_t len>
template <size_t len2>
fix_int<len>& fix_int<len>::operator=(fix_int<len2>&& x) {
  sign = x.sign;
  if (len <= len2) {
    std::memcpy(data, x.data, sizeof(data));
  } else {
    std::memcpy(data, x.data, x.data_len * 8);
    std::memset(data + x.data_len, 0x00, 8 * (data_len - x.data_len));
  }
  return *this;
}

template <size_t len>
char fix_int<len>::to_char(int x) {
  if (x < 10) [[likely]] {
    return x + '0';
  } else {
    return x - 10 + 'a';
  }
}

template <size_t len>
int fix_int<len>::from_char(char ch) {
  if (ch <= '9') [[likely]] {  // '9' < 'a'.
    return ch - '0';
  } else {
    return ch - 'a' + 10;
  }
}

template <size_t len>
std::string fix_int<len>::to_hex_string() {
  std::string ret;
  ret.resize(len / 4);
  for (int i = 0; i < data_len; ++i) {
    uint64 tmp = data[data_len - 1 - i];
    int j = 16 * i + 15;
    for (int k = 0; k <= 15; ++k) {
      ret[j] = to_char(tmp & 0xF);
      tmp >>= 4;
      j--;
    }
  }
  int i = 0;
  while (ret[i] == '0' && i + 1 < len / 4) {
    i++;
  }
  ret = ret.substr(i);
  if (sign == true) {
    ret = "-" + ret;
  }
  return ret;
}

template <size_t len>
void fix_int<len>::from_hex_string(std::string str) {
  sign = false;
  if (str[0] == '-') {
    sign = true;
    str = str.substr(1);
  }
  int length = str.length();
  int chunks = (length - 1) / 16;
  int now = chunks;  // now position of data
  int start_size = length - chunks * 16;
  uint64 tmp = 0;
  for (int i = 0; i < start_size; ++i) {
    tmp = (tmp << 4) | from_char(str[i]);
  }
  data[now] = tmp;
  now--;
  for (int i = start_size; i < length; i += 16) {
    uint64 tmp = 0;
    for (int j = 0; j <= 15; ++j) {
      tmp = (tmp << 4) | from_char(str[i + j]);
    }
    data[now] = tmp;
    now--;
  }
  if (is_zero() && sign) {
    sign = false;
  }
}

template <size_t len>
inline fix_int<len> operator+(fix_int<len> const& x) noexcept {
  return x;
}
template <size_t len>
inline fix_int<len> operator-(fix_int<len> x) {
  x.sign ^= true;
  return x;
}

template <size_t len>
inline fix_int<len> operator++(fix_int<len>& x) {
  return x += 1;
}
template <size_t len>
inline fix_int<len> operator++(fix_int<len>& x, int t) {
  fix_int<len> ret = x;
  x += 1;
  return ret;
}
template <size_t len>
inline fix_int<len> operator--(fix_int<len>& x) {
  return x -= 1;
}
template <size_t len>
inline fix_int<len> operator--(fix_int<len>& x, int t) {
  fix_int<len> ret = x;
  x -= 1;
  return ret;
}

template <size_t len1, size_t len2>
inline fix_int<len1>& operator+=(fix_int<len1>& lhs, fix_int<len2> const& rhs) {
  return lhs = lhs + rhs;
}
template <size_t len1, size_t len2>
inline fix_int<len1>& operator-=(fix_int<len1>& lhs, fix_int<len2> const& rhs) {
  return lhs = lhs - rhs;
}
template <size_t len1, size_t len2>
inline fix_int<len1>& operator*=(fix_int<len1>& lhs, fix_int<len2> const& rhs) {
  return lhs = lhs * rhs;
}

template <size_t len1, size_t len2>
inline fix_int<len1>& operator/=(fix_int<len1>& lhs, fix_int<len2> const& rhs) {
  return lhs = lhs / rhs;
}
template <size_t len1, std::unsigned_integral T>
inline fix_int<len1>& operator/=(fix_int<len1>& lhs, T rhs) {
  return lhs = lhs / rhs;
}
template <size_t len1, std::signed_integral T>
inline fix_int<len1>& operator/=(fix_int<len1>& lhs, T rhs) {
  return lhs = lhs / rhs;
}
template <size_t len1, size_t len2>
inline fix_int<len1>& operator%=(fix_int<len1>& lhs, fix_int<len2> const& rhs) {
  return lhs = lhs % rhs;
}
template <size_t len1, std::unsigned_integral T>
inline fix_int<len1>& operator%=(fix_int<len1>& lhs, T rhs) {
  return lhs = lhs % rhs;
}
template <size_t len1, std::signed_integral T>
inline fix_int<len1>& operator%=(fix_int<len1>& lhs, T rhs) {
  return lhs = lhs % rhs;
}

// template <size_t len>
// inline fix_int<len>::operator bool() const {
//   return !is_zero();
// }
// template <size_t len>
// inline fix_int<len>::operator uint64() const {
//   return data[0];
// }

template <size_t len1, size_t len2>
fix_int<std::max(len1, len2)> operator&(fix_int<len1> const& lhs,
                                        fix_int<len2> const& rhs) {
  fix_int<std::max(len1, len2)> ret = lhs;
  for (int i = 0; i < rhs.data_len; ++i) {
    ret.data[i] &= rhs.data[i];
  }
  return ret;
}
template <size_t len1, size_t len2>
fix_int<std::max(len1, len2)> operator|(fix_int<len1> const& lhs,
                                        fix_int<len2> const& rhs) {
  fix_int<std::max(len1, len2)> ret = lhs;
  for (int i = 0; i < rhs.data_len; ++i) {
    ret.data[i] |= rhs.data[i];
  }
  return ret;
}
template <size_t len1, size_t len2>
fix_int<std::max(len1, len2)> operator^(fix_int<len1> const& lhs,
                                        fix_int<len2> const& rhs) {
  fix_int<std::max(len1, len2)> ret = lhs;
  for (int i = 0; i < rhs.data_len; ++i) {
    ret.data[i] ^= rhs.data[i];
  }
  return ret;
}
template <size_t len>
fix_int<len> operator~(fix_int<len> x) {
  for (int i = 0; i < x.data_len; ++i) {
    x.data[i] = ~x.data[i];
  }
  return x;
}

template <size_t len1, size_t len2>
inline fix_int<len1>& operator&=(fix_int<len1>& lhs, fix_int<len2> const& rhs) {
  return lhs = lhs & rhs;
}
template <size_t len1, size_t len2>
inline fix_int<len1>& operator|=(fix_int<len1>& lhs, fix_int<len2> const& rhs) {
  return lhs = lhs | rhs;
}
template <size_t len1, size_t len2>
inline fix_int<len1>& operator^=(fix_int<len1>& lhs, fix_int<len2> const& rhs) {
  return lhs = lhs ^ rhs;
}

template <size_t len1>
fix_int<len1> operator>>(fix_int<len1> const& x, size_t digits) {
  size_t chunks = digits / 64, offset = digits % 64;
  fix_int<len1> ret;
  ret.sign = x.sign;
  for (int i = x.data_len - 1 - chunks; i >= 0; --i) {
    ret.data[i] = x.data[i + chunks] >> offset;
    if (offset != 0 && i < x.data_len - 1 - chunks) {
      ret.data[i] |= x.data[i + chunks + 1] << (64 - offset);
    }
  }
  return ret;
}
template <size_t len1>
fix_int<len1> operator<<(fix_int<len1> const& x, size_t digits) {
  size_t chunks = digits / 64, offset = digits % 64;
  fix_int<len1> ret;
  ret.sign = x.sign;
  for (int i = chunks; i < x.data_len; ++i) {
    ret.data[i] = x.data[i - chunks] << offset;
    if (offset != 0 && i > chunks) {
      ret.data[i] |= x.data[i - chunks - 1] >> (64 - offset);
    }
  }
  return ret;
}

template <size_t len1>
inline fix_int<len1>& operator<<=(fix_int<len1>& lhs, size_t digits) {
  return lhs = lhs << digits;
}
template <size_t len1>
inline fix_int<len1>& operator>>=(fix_int<len1>& lhs, size_t digits) {
  return lhs = lhs >> digits;
}

template <size_t len>
template <size_t len1>
void fix_int<len>::add_data(fix_int<len1> const& rhs) {
  static_assert(len >= len1,
                "fix_int::add_data: len should be no less than len1.");
  uint64 carry = 0;
  for (int i = 0; i < rhs.data_len; ++i) {
    if (carry) {
      data[i]++;
      if (data[i] == 0) [[unlikely]] {
        carry = 1;
      } else {
        carry = 0;
      }
    }
    data[i] += rhs.data[i];
    if (data[i] < rhs.data[i]) {
      carry = 1;
    }
  }
  if (carry) [[unlikely]] {
    int i = rhs.data_len;
    while (true) {
      data[i]++;
      if (data[i] != 0) [[likely]] {
        break;
      }
      i++;
    }
  }
}

template <size_t len>
template <size_t len1>
// clang-format off
requires(len >= len1) 
int fix_int<len>::compare(fix_int<len1> const& rhs) const {
  // clang-format on
  for (int i = data_len - 1; i >= rhs.data_len; --i) {
    if (data[i] != 0) {
      return 1;
    }
  }
  for (int i = rhs.data_len - 1; i >= 0; --i) {
    if (data[i] == rhs.data[i]) {
      continue;
    }
    return data[i] > rhs.data[i] ? 1 : -1;
  }
  return 0;
}
template <size_t len>
template <size_t len1>
// clang-format off
// this 2 > 1 is for coloring for template stuff.
requires(len < len1 && 2 > 1)
inline int fix_int<len>::compare(fix_int<len1> const& rhs) const {
  // clang-format on
  return -rhs.compare(*this);
}

template <size_t len>
template <size_t len1>
void fix_int<len>::subtract_data(fix_int<len1> const& rhs) {
  uint64 carry = 0;
  for (int i = 0; i < data_len; ++i) {
    if (carry) {
      carry = 0;
      if (data[i] == 0) [[unlikely]] {
        carry = 1;
      }
      data[i]--;
    }
    if (i >= rhs.data_len) {
      continue;
    }
    uint64 tmp = data[i];
    data[i] -= rhs.data[i];
    if (data[i] > tmp) {
      carry = 1;
    }
  }
}

template <size_t len, size_t len1>
fix_int<std::max(len, len1) + 64> operator+(fix_int<len> const& lhs,
                                            fix_int<len1> const& rhs) {
  fix_int<std::max(len, len1) + 64> ret;
  if (lhs.sign == rhs.sign) {
    ret = lhs;
    ret.add_data(rhs);
  } else {
    int comp = lhs.compare(rhs);
    if (comp > 0) {
      ret = lhs;
      ret.subtract_data(rhs);
    } else if (comp < 0) {
      ret = rhs;
      ret.subtract_data(lhs);
    }
  }
  return ret;
}

template <size_t len, size_t len1>
fix_int<std::max(len, len1) + 64> operator-(fix_int<len> const& lhs,
                                            fix_int<len1> const& rhs) {
  fix_int<std::max(len, len1) + 64> ret;
  if (lhs.sign != rhs.sign) {
    ret = lhs;
    ret.add_data(rhs);
  } else {
    int cmp = lhs.compare(rhs);
    if (cmp > 0) {
      ret = lhs;
      ret.subtract_data(rhs);
    } else if (cmp < 0) {
      ret = rhs;
      ret.subtract_data(lhs);
      ret.sign ^= true;
    }
  }
  return ret;
}

template <size_t len>
template <size_t len1, size_t len2>
void fix_int<len>::long_multiplication(fix_int<len1> const& lhs,
                                       fix_int<len2> const& rhs) {
  uint128 tmp[len / 64];
  uint128 carry[len / 64];
  std::memset(tmp, 0x00, sizeof(tmp));
  std::memset(carry, 0x00, sizeof(carry));
  int n = lhs.first_significant_chunk(), m = rhs.first_significant_chunk();
  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <= m; ++j) {
      uint128 product = static_cast<uint128>(lhs.data[i]) * rhs.data[j];
      tmp[i + j] += product;
      if (tmp[i + j] < product) {
        carry[i + j]++;
      }
    }
  }
  for (int i = 1; i <= n + m + 1; ++i) {
    if (i >= 2) {
      tmp[i] += carry[i - 2];
      if (tmp[i] < carry[i - 2]) {
        carry[i]++;
      }
    }
    tmp[i] += tmp[i - 1] >> 64;
    if (tmp[i] < tmp[i - 1] >> 64) {
      carry[i]++;
    }
  }
  for (int i = 0; i <= n + m + 1; ++i) {
    data[i] = tmp[i];
  }
}
template <size_t len>
template <size_t len1, size_t len2>
void fix_int<len>::karatsuba(fix_int<len1> const& lhs,
                             fix_int<len2> const& rhs) {
  static_assert(len1 >= 512 && len2 >= 512,
                "fix_int::karatuba needs len1, len2 >= 512.");
  constexpr size_t bits = (std::min(len1, len2) / 2) / 64 * 64;
  constexpr size_t size = std::max(len1, len2) - bits;
  fix_int<size> a, b, c, d;
  std::memcpy(a.data, lhs.data + bits / 64, (len1 - bits) / 8);
  std::memcpy(b.data, lhs.data, bits / 8);
  std::memcpy(c.data, rhs.data + bits / 64, (len2 - bits) / 8);
  std::memcpy(d.data, rhs.data, bits / 8);
  // clang-format off
  fix_int<size * 2 + 64> tmp1 = a * c, tmp2 = (a + b) * (c + d), tmp3 = b * d;
  // clang-format on
  tmp2 -= tmp1;
  tmp2 -= tmp3;
  std::memcpy(data + bits / 64, tmp1.data, size * 2 / 8);
  add_data(tmp2);
  *this <<= bits;
  add_data(tmp3);
}

/**
 * @brief The implementation is based on Marco Bodrato, Torwards Optimal 
 * Toom-Cook Multiplication for Univariate and Multivariate Polynomials in
 * Characteristic 2 and 0, Arithmetic of Finite Fields, 2007.
 */
template <size_t len>
template <size_t len1, size_t len2>
void fix_int<len>::toom_cook(fix_int<len1> const& lhs,
                             fix_int<len2> const& rhs) {
  static_assert(len1 >= 1024 && len2 >= 1024,
                "fix_int::karatuba needs len1, len2 >= 1024.");
  constexpr size_t bits = (std::min(len1, len2) / 3) / 64 * 64;
  constexpr size_t size = std::max(len1, len2) - 2 * bits;
  fix_int<size> u2, u1, u0, v2, v1, v0;
  std::memcpy(u2.data, lhs.data + bits / 64 * 2, (len1 - bits * 2) / 8);
  std::memcpy(u1.data, lhs.data + bits / 64, bits / 8);
  std::memcpy(u0.data, lhs.data, bits / 8);
  std::memcpy(v2.data, rhs.data + bits / 64 * 2, (len2 - bits * 2) / 8);
  std::memcpy(v1.data, rhs.data + bits / 64, bits / 8);
  std::memcpy(v0.data, rhs.data, bits / 8);
  // The upper bound is calculated carefully by interval arthimetic.
  // Note that the result will be in an order of w4, w2, w1, w3, and w0.
  fix_int<2 * size + 64> w0, w1, w2, w3, w4;
  fix_int<size + 64> w0_, w1_, w2_, w3_, w4_;
  w0_ = u2 + u0;
  w4_ = v2 + v0;
  w2_ = w0_ - u1;
  w1_ = w4_ - v1;
  w0_ = w0_ + u1;
  w4_ = w4_ + v1;
  w3 = w2_ * w1_;
  w1 = w0_ * w4_;
  // TODO: Compare the efficiency of *2 and <<1.
  w0_ = (w0_ + u2) * 2 - u0;
  w4_ = (w4_ + v2) * 2 - v0;
  w2 = w0_ * w4_;
  w0 = u0 * v0;
  w4 = u2 * v2;
  w2 = (w2 - w3) / 3;
  w3 = (w1 - w3) / 2;
  w1 = w1 - w0;
  w2 = (w2 - w1) / 2 - w4 * 2;
  w1 = w1 - w3 - w4;
  w3 = w3 - w2;
  // It can be shown that w0, w4 >= 0.
  std::memcpy(data + bits / 64, w4.data, size * 2 / 8);
  *this += w2;
  *this <<= bits;
  *this += w1;
  *this <<= bits;
  *this += w3;
  *this <<= bits;
  add_data(w0);
}

/**
 * @brief The implementation is based on 
 * https://en.algorithmica.org/hpc/number-theory/montgomery/ and 
 * https://codeforces.com/blog/entry/129600 .
 */
template <size_t len>
inline uint64 fix_int<len>::number_theoretic_transform_t::montgomery_t::reduce(
    uint128 const& x) {
  uint64 q = x * inv_mod;
  uint64 m = (static_cast<uint128>(q) * mod) >> 64;
  uint64 res = (x >> 64) + mod - m;
  return res >= mod ? res - mod : res;
}
template <size_t len>
inline uint64
fix_int<len>::number_theoretic_transform_t::montgomery_t::multiply(uint64 x,
                                                                   uint64 y) {
  return reduce(static_cast<uint128>(x) * y);
  // return static_cast<uint128>(x) * y % mod;
}
template <size_t len>
inline uint64
fix_int<len>::number_theoretic_transform_t::montgomery_t::to_montgomery(
    uint64 x) {
  return (static_cast<uint128>(x) << 64) % mod;
  // return x;
}
template <size_t len>
inline uint64
fix_int<len>::number_theoretic_transform_t::montgomery_t::from_montgomery(
    uint64 x) {
  return reduce(x);
  // return x;
}

template <size_t len>
std::vector<uint64> fix_int<len>::number_theoretic_transform_t::calc_bit_inv() {
  std::vector<uint64> ret;
  ret.resize(1 << 15);
  ret[1] = 1 << 14;
  for (int i = 2; i < (1 << 15); ++i) {
    ret[i] = (ret[i >> 1] >> 1) | ((i & 1) << 14);
  }
  return ret;
}

template <size_t len>
inline uint64 fix_int<len>::number_theoretic_transform_t::bit_inverse(
    uint64 x) {
  return bit_inv[x >> 15] | (bit_inv[x & 0x7fff] << 15);
}

template <size_t len>
std::vector<uint64> fix_int<len>::number_theoretic_transform_t::bit_inv =
    std::vector<uint64>(1 << 15, 0);

template <size_t len>
uint64 fix_int<len>::number_theoretic_transform_t::modular_pow(uint64 x,
                                                               uint64 y) {
  uint128 res = 1;
  uint128 now = x;
  while (y) {
    if (y & 1) {
      res = res * now % mod;
    }
    now = now * now % mod;
    y >>= 1;
  }
  return res;
}

template <size_t len>
inline uint64 fix_int<len>::number_theoretic_transform_t::modular_inv(
    uint64 x) {
  return modular_pow(x, mod - 2);
}

template <size_t len>
void fix_int<len>::number_theoretic_transform_t::transform(
    std::vector<uint64>& a, bool on) {
  int n = a.size();
  int u = 1, v = 0;
  while (u != n) {
    u *= 2;
    v++;
  }
  for (int i = 0; i < n; ++i) {
    int rev_i = bit_inverse(i) >> (30 - v);
    if (i < rev_i) {
      std::swap(a[i], a[rev_i]);
    }
  }
  for (int h = 2; h <= n; h <<= 1) {
    uint64 omega_n = modular_pow(proot, (mod - 1) / h);
    if (on) {
      omega_n = modular_inv(omega_n);
    }
    omega_n = montgomery_t::to_montgomery(omega_n);
    for (int j = 0; j < n; j += h) {
      uint64 omega = montgomery_t::to_montgomery(1);
      for (int k = j; k < j + h / 2; ++k) {
        uint64 u = a[k];
        uint64 v = montgomery_t::multiply(omega, a[k + h / 2]);
        a[k] = u + v;
        a[k + h / 2] = u + mod - v;
        if (a[k] >= mod) {
          a[k] -= mod;
        }
        if (a[k + h / 2] >= mod) {
          a[k + h / 2] -= mod;
        }
        omega = montgomery_t::multiply(omega, omega_n);
      }
    }
  }
  if (on) {
    uint64 inv_n = montgomery_t::to_montgomery(modular_inv(n));
    for (int i = 0; i < n; ++i) {
      a[i] = montgomery_t::multiply(a[i], inv_n);
    }
  }
}

template <size_t len>
template <size_t len1, size_t len2>
void fix_int<len>::number_theoretic_transform(fix_int<len1> const& lhs,
                                              fix_int<len2> const& rhs) {
  if (number_theoretic_transform_t::bit_inv[1] == 0) {
    number_theoretic_transform_t::bit_inv =
        number_theoretic_transform_t::calc_bit_inv();
  }
  size_t n = 64;
  while (n < len1 + len2) {
    n *= 2;
  }
  n /= 16;
  std::vector<uint64> a, b;
  a.resize(n);
  b.resize(n);
  for (int i = 0; i < lhs.data_len; ++i) {
    a[4 * i] = lhs.data[i] & 0xffff;
    a[4 * i + 1] = (lhs.data[i] >> 16) & 0xffff;
    a[4 * i + 2] = (lhs.data[i] >> 32) & 0xffff;
    a[4 * i + 3] = (lhs.data[i] >> 48) & 0xffff;
  }
  for (int i = 0; i < rhs.data_len; ++i) {
    b[4 * i] = rhs.data[i] & 0xffff;
    b[4 * i + 1] = (rhs.data[i] >> 16) & 0xffff;
    b[4 * i + 2] = (rhs.data[i] >> 32) & 0xffff;
    b[4 * i + 3] = (rhs.data[i] >> 48) & 0xffff;
  }
  for (int i = 0; i < n; ++i) {
    a[i] = number_theoretic_transform_t::montgomery_t::to_montgomery(a[i]);
    b[i] = number_theoretic_transform_t::montgomery_t::to_montgomery(b[i]);
  }
  number_theoretic_transform_t::transform(a, false);
  number_theoretic_transform_t::transform(b, false);
  for (int i = 0; i < n; ++i) {
    a[i] = number_theoretic_transform_t::montgomery_t::multiply(a[i], b[i]);
  }
  number_theoretic_transform_t::transform(a, true);
  for (int i = 0; i < n; ++i) {
    a[i] = number_theoretic_transform_t::montgomery_t::from_montgomery(a[i]);
  }
  uint64 carry = 0;
  for (int i = 0; i < data_len && 4 * i < n; ++i) {
    data[i] = carry;
    carry = 0;
    for (int j = 0; j < 4; ++j) {
      uint64 tmp = a[4 * i + j] << (16 * j);
      data[i] += tmp;
      if (data[i] < tmp) {
        carry++;
      }
    }
    if (i > 0) {
      for (int j = 1; j <= 3; ++j) {
        uint64 tmp = a[4 * i - j] >> (16 * j);
        data[i] += tmp;
        if (data[i] < tmp) {
          carry++;
        }
      }
    }
  }
  if (n / 4 < data_len) {
    data[n / 4] = carry;
    data[n / 4] += a[n - 1] >> 16;
    data[n / 4] += a[n - 2] >> 32;
    data[n / 4] += a[n - 3] >> 48;
  }
}

template <size_t len>
bool fix_int<len>::is_zero() const {
  for (int i = 0; i < len / 64; ++i) {
    if (data[i] != 0) {
      return false;
    }
  }
  return true;
}

template <size_t len1, size_t len2>
inline fix_int<len1 + len2> operator*(fix_int<len1> const& lhs,
                                      fix_int<len2> const& rhs) {
  fix_int<len1 + len2> ret;
  if constexpr (true || (len1 >= 2048 && len2 >= 2048)) {
    ret.number_theoretic_transform(lhs, rhs);
  } else if constexpr (len1 >= 1024 && len2 >= 1024) {
    ret.toom_cook(lhs, rhs);
  } else if constexpr (len1 >= 512 && len2 >= 512) {
    ret.karatsuba(lhs, rhs);
  } else {
    ret.long_multiplication(lhs, rhs);
  }
  ret.sign = !ret.is_zero() && (lhs.sign ^ rhs.sign);
  return ret;
}
template <size_t len1>
inline fix_int<len1 + 64> operator*(fix_int<len1> lhs, uint64 rhs) {
  return lhs * fix_int<64>(rhs);
}

template <size_t len1, std::unsigned_integral T>
std::pair<fix_int<len1>, fix_int<64>> division(fix_int<len1> const& lhs,
                                               T rhs) {
  if (rhs == 0) [[unlikely]] {
    throw std::invalid_argument("yuki::division: divided by 0.");
  }
  uint128 now = 0;
  fix_int<len1> ret;
  for (int i = len1 / 64 - 1; i >= 0; --i) {
    now = now << 64 | lhs.data[i];
    uint128 quotient = now / rhs;
    ret.data[i] = quotient;
    now = now - quotient * rhs;
  }
  fix_int<64> ret2 = static_cast<uint64>(now);
  if (lhs.sign && !ret.is_zero()) {
    ret.sign = true;
  }
  if (now != 0) {
    ret2.sign = lhs.sign;
  }
  return {ret, ret2};
}

template <size_t len1, std::unsigned_integral T>
inline fix_int<len1> operator/(fix_int<len1> const& lhs, T rhs) {
  return division(lhs, rhs).first;
}

template <size_t len1, std::unsigned_integral T>
fix_int<64> operator%(fix_int<len1> const& lhs, T rhs) {
  if (rhs == 0) [[unlikely]] {
    throw std::invalid_argument("yuki::division: divided by 0.");
  }
  uint128 now = 0;
  for (int i = len1 / 64 - 1; i >= 0; --i) {
    now = (now << 64 | lhs.data[i]) % rhs;
  }
  fix_int<64> ret2 = static_cast<uint64>(now);
  if (now != 0) [[likely]] {
    ret2.sign = lhs.sign;
  }
  return ret2;
}

template <size_t len1, std::signed_integral T>
std::pair<fix_int<len1>, fix_int<64>> division(fix_int<len1> const& lhs,
                                               T rhs) {
  if (rhs == 0) {
    throw std::invalid_argument("yuki::division: divided by 0.");
  }
  uint64 rhs2 = std::abs(rhs);
  auto [ret, now] = division(lhs, rhs2);
  if (lhs.sign && !now.is_zero()) {
    now.sign ^= true;
  }
  return {ret, now};
}

template <size_t len1, std::signed_integral T>
inline fix_int<len1> operator/(fix_int<len1> const& lhs, T rhs) {
  return division(lhs, rhs).first;
}

template <size_t len1, std::signed_integral T>
inline fix_int<64> operator%(fix_int<len1> const& lhs, T rhs2) {
  if (rhs2 == 0) {
    throw std::invalid_argument("yuki::division: divided by 0.");
  }
  uint64 rhs = rhs2;
  uint128 now = 0;
  for (int i = len1 / 64 - 1; i >= 0; --i) {
    now = (now << 64 | lhs.data[i]) % rhs;
  }
  fix_int<64> ret2 = static_cast<uint64>(now);
  if (now != 0) {
    ret2.sign = lhs.sign;
  }
  return ret2;
}

template <size_t len>
size_t fix_int<len>::first_significant_chunk() const {
  for (int i = data_len - 1; i >= 0; --i) {
    if (data[i]) {
      return i;
    }
  }
  return -1;
}

template <size_t len>
template <size_t len1, size_t len2>
void fix_int<len>::schoolbook_division(fix_int<len1> lhs, fix_int<len2> rhs) {
  if (rhs.is_zero()) [[unlikely]] {
    throw std::invalid_argument("yuki::division: divided by 0.");
  }
  lhs.sign = rhs.sign = false;
  int n = lhs.first_significant_chunk(), m = rhs.first_significant_chunk();
  if (m == 0) {
    *this = lhs / rhs.data[0];
    return;
  }
  constexpr uint128 B = static_cast<uint128>(1) << 64;
  if (rhs.data[m] < B / 2) {
    lhs *= fix_int<64>(static_cast<uint64>(B / (rhs.data[m] + 1)));
    rhs *= fix_int<64>(static_cast<uint64>(B / (rhs.data[m] + 1)));
    n = lhs.first_significant_chunk();
  }
  fix_int<len1> aligned_rhs = rhs << ((n - m) * 64);
  for (int i = n; i >= m; --i) {
    uint128 lhs_approx = (i == n ? 0 : lhs.data[i + 1]) * B + lhs.data[i];
    uint128 q = std::min(lhs_approx / rhs.data[m], B - 1);
    if (q == 0) {
      aligned_rhs >>= 64;
      continue;
    }
    uint128 t = lhs_approx - q * rhs.data[m];
    uint128 tmp2 = static_cast<uint128>(q) * rhs.data[m - 1];
    if (tmp2 > lhs.data[i - 1] && B * t < tmp2 - lhs.data[i - 1]) {
      q--;
    }
    fix_int<len1> tmp = aligned_rhs * static_cast<uint64>(q);
    if (tmp.compare(lhs) > 0) {
      q--;
      tmp -= aligned_rhs;
    }
    lhs -= tmp;
    data[i - m] = q;
    aligned_rhs >>= 64;
  }
}

template <size_t len>
inline constexpr int fix_int<len>::newton_raphson_times() {
  long double p = len + 1;
  p /= std::log2l(17);
  return std::ceil(std::log2l(p));
}

template <size_t len>
fix_int<len>::fix_point_t::fix_point_t() {
  value = fix_int<len>();
  offset = 0;
}

template <size_t len>
template <size_t len1, size_t len2>
typename fix_int<len1 + len2>::fix_point_t fix_int<len>::add_fix_point(
    typename fix_int<len1>::fix_point_t const& lhs,
    typename fix_int<len2>::fix_point_t const& rhs) {
  typename fix_int<len1 + len2>::fix_point_t ret;
  if (lhs.offset < rhs.offset) {
    ret.value = lhs;
    ret.value <<= (rhs.offset - lhs.offset);
    ret.value += rhs.value;
    ret.offset = rhs.offset;
  } else {
    ret.value = rhs;
    ret.value <<= (lhs.offset - rhs.offset);
    ret.value += lhs.value;
    ret.offset = lhs.offset;
  }
  return ret;
}

template <size_t len>
template <size_t len1, size_t len2>
inline typename fix_int<len1 + len2>::fix_point_t fix_int<len>::mul_fix_point(
    typename fix_int<len1>::fix_point_t const& lhs,
    typename fix_int<len2>::fix_point_t const& rhs) {
  typename fix_int<len1 + len2>::fix_point_t ret;
  ret.value = lhs.value * rhs.value;
  ret.offset = lhs.offset + rhs.offset;
  return ret;
}

template <size_t len1>
typename fix_int<len1>::fix_point_t trunc_fix_point(
    typename fix_int<len1>::fix_point_t const& x, int off) {
  assert(off <= x.offset);
  int start = x.offset - off;
  x.value.data[start / 64] >>= start % 64;
  x.value.data[start / 64] <<= start % 64;
  std::memset(x.value.data, 0x00, start / 64 * 8);
}

template <size_t len>
inline bool fix_int<len>::get_bit(int index) {
  return (data[index / 64] >> (index % 64)) & 1;
}

// template <size_t len>
// template <size_t len1, size_t len2>
// void fix_int<len>::newton_raphson_division(fix_int<len1> const& lhs,
//                                            fix_int<len2> const& rhs) {
//   if (rhs.is_zero()) [[unlikely]] {
//     throw std::invalid_argument("yuki::newton_raphson_division: divided by 0.");
//   }
//   int n1 = lhs.first_significant_chunk(), n2 = rhs.first_significant_chunk();
//   int m1 = std::countl_zero(lhs.data[n1]) + n1 * 64;
//   int m2 = std::countl_zero(lhs.data[n2]) + n2 * 64;
//   if (m1 < m2) {
//     *this = 0;
//     return;
//   }
//   if (m2 <= 5) {
//     *this = lhs / rhs.data[0];
//     return;
//   }
//   int digits = m1 + 1 - m2;
//   int size = (2 * digits + 63) / 64 * 64;
//   typename fix_int<(len1 + len2) * 2 + 64>::fix_point_t z, s, t, u, w;
//   z.value = 32 / (4 * rhs.get_bit(m2) + 2 * rhs.get_bit(m2 - 1) +
//                   rhs.get_bit(m2 - 2));
//   z.offset = 2;
//   for (int k = 1; k < digits; k *= 2) {
//     s = mul_fix_point(z, z);
//     t.value = rhs;
//     t.value.data = false;
//     t.offset = m2;
//     trunc_fix_point(t, 2 * k + 3);
//     u = mul_fix_point(t, s);
//     trunc_fix_point(u, 2 * k + 1);
//     w = add_fix_point(z, z);
//     u.data.sign ^= true;
//     z = add_fix_point(w, u);
//   }
//   auto tmp = lhs * z.value;
//   tmp >>= z.offset + m2;
//   *this = tmp;
//   this->sign = false;
// }

template <size_t len1, size_t len2>
std::pair<fix_int<len1>, fix_int<std::min(len1, len2)>> division(
    fix_int<len1> const& lhs, fix_int<len2> const& rhs) {
  fix_int<len1> quotient;
  // quotient.newton_raphson_division(lhs, rhs);
  quotient.schoolbook_division(lhs, rhs);
  if (!quotient.is_zero()) {
    quotient.sign = lhs.sign ^ rhs.sign;
  }
  return {quotient, lhs - rhs * quotient};
}

template <size_t len1, size_t len2>
fix_int<len1> operator/(fix_int<len1> const& lhs, fix_int<len2> const& rhs) {
  fix_int<len1> quotient;
  // quotient.newton_raphson_division(lhs, rhs);
  quotient.schoolbook_division(lhs, rhs);
  if (!quotient.is_zero()) {
    quotient.sign = lhs.sign ^ rhs.sign;
  }
  return quotient;
}

template <size_t len1, size_t len2>
fix_int<std::min(len1, len2)> operator%(fix_int<len1> const& lhs,
                                        fix_int<len2> const& rhs) {
  return division(lhs, rhs).second;
}

}  // namespace yuki

#endif  // #ifndef YUKI_LIB_FIX_INT_H_
