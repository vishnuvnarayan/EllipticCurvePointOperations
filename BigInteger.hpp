/**
    BigInteger.hpp
    BigInteger class (header file) to manipulate arbitrary-length signed
integers.
    Author: Vishnu V Narayan
**/

#pragma once

#include <climits>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <utility>
#include <vector>

const int32_t BI_BASE = 1000000;
const int32_t BI_BASE_DIGITS = 6;

// forward declarations:
class BigInteger;
bool BI_Miller_Rabin(BigInteger p, int32_t it = 32);

class BigInteger {

  // Stores arbitrarily large integer absolute value
  std::vector<int32_t> a;
  int32_t sign;

public:
  // Helper functions
  void read(const std::string &);
  void trim();
  friend std::vector<int64_t> karatsuba(const std::vector<int64_t> &,
                                        const std::vector<int64_t> &);
  friend std::vector<int32_t> BI_Helper_DecBin(const BigInteger &);
  friend BigInteger BI_Helper_BinDec(const std::vector<int32_t> &);
  friend std::vector<int32_t> BI_Helper_Subtract(const std::vector<int32_t> &,
                                                 const std::vector<int32_t> &);
  friend bool BI_Helper_Less(std::vector<int32_t>, std::vector<int32_t>);
  friend int32_t BI_Helper_CompareAbs(const BigInteger &, const BigInteger &);
  friend std::pair<BigInteger, BigInteger> BI_Helper_Divide(const BigInteger &,
                                                            const BigInteger &);

  // Constructors
  BigInteger();
  BigInteger(int32_t);
  BigInteger(int64_t);
  BigInteger(const BigInteger &);
  BigInteger(const std::string &);

  // operator=
  BigInteger &operator=(int32_t);
  BigInteger &operator=(int64_t);
  BigInteger &operator=(const BigInteger &);

  // istream and ostream operators
  friend std::istream &operator>>(std::istream &, BigInteger &);
  friend std::ostream &operator<<(std::ostream &, const BigInteger &);

  // Relational operators
  bool operator<(const BigInteger &) const;
  bool operator>(const BigInteger &) const;
  bool operator<=(const BigInteger &) const;
  bool operator>=(const BigInteger &) const;
  bool operator==(const BigInteger &) const;
  bool operator!=(const BigInteger &) const;

  // int64_t relational operators
  bool operator<(const int64_t) const;
  bool operator>(const int64_t) const;
  bool operator<=(const int64_t) const;
  bool operator>=(const int64_t) const;
  bool operator==(const int64_t) const;
  bool operator!=(const int64_t) const;

  // Compound arithmetic operators
  BigInteger &operator+=(const BigInteger &);
  BigInteger &operator-=(const BigInteger &);
  BigInteger &operator*=(const BigInteger &);
  BigInteger &operator/=(const BigInteger &);
  BigInteger &operator%=(const BigInteger &);

  // Arithmetic operators (non-members)
  friend const BigInteger operator+(BigInteger, const BigInteger &);
  friend const BigInteger operator-(BigInteger, const BigInteger &);
  friend const BigInteger operator*(BigInteger, const BigInteger &);
  friend const BigInteger operator/(BigInteger, const BigInteger &);
  friend const BigInteger operator%(BigInteger, const BigInteger &);

  // int64_t compound arithmetic operators
  BigInteger &operator+=(int64_t);
  BigInteger &operator-=(int64_t);
  BigInteger &operator*=(int64_t);
  BigInteger &operator/=(int64_t);
  BigInteger &operator%=(int64_t);

  // int64_t arithmetic operators (non-members)
  friend const BigInteger operator+(BigInteger, int64_t);
  friend const BigInteger operator-(BigInteger, int64_t);
  friend const BigInteger operator*(BigInteger, int64_t);
  friend const BigInteger operator/(BigInteger, int64_t);
  friend const BigInteger operator%(BigInteger, int64_t);

  // pre- and post-increment and decrement operators
  BigInteger &operator++();
  BigInteger operator++(int32_t);
  BigInteger &operator--();
  BigInteger operator--(int32_t);

  // Miscellaneous
  int64_t to_int64() const;
  bool is_even() const;
  bool is_zero() const;
  bool is_negative() const;
  friend BigInteger abs(BigInteger &);
  friend BigInteger gcd(BigInteger, BigInteger);
  friend BigInteger lcm(BigInteger, BigInteger);
  friend BigInteger BI_FastExp(BigInteger, BigInteger, BigInteger);
  friend BigInteger BI_FastExp(BigInteger, int64_t, BigInteger);
  friend int32_t BI_FastExp(BigInteger, int64_t, int32_t);
  friend BigInteger BI_ModInv(BigInteger, BigInteger);
  friend BigInteger BI_Residue(BigInteger, BigInteger);
  friend bool BI_Miller_Rabin(BigInteger, int32_t);
};
