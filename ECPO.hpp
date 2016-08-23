#pragma once

#include "BigInteger.hpp"
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

class EC_Point {
public:
  BigInteger x;
  BigInteger y;
  bool zero;

  // Constructors
  EC_Point();
  EC_Point(int32_t);
  EC_Point(int64_t);
  EC_Point(const BigInteger &, const BigInteger &);
  EC_Point(const EC_Point &);
  EC_Point(std::pair<BigInteger, BigInteger> &);

  // ostream operator
  friend std::ostream &operator<<(std::ostream &, const EC_Point &);
};

class EC_Curve {
private:
  bool valid;
  bool warning;

public:
  BigInteger p;
  BigInteger A;
  BigInteger B;

  // Constructors
  EC_Curve();
  EC_Curve(const BigInteger &, const BigInteger &, const BigInteger &);

  // EC_Validate - validate the curve for point operations (p prime,
  // discriminant not 0)
  void EC_Validate();

  // EC_Validate_Point - returns true if the parameter point is on the curve
  bool EC_Validate_Point(EC_Point);

  // EC_Add - add two points on the curve
  EC_Point EC_Add(EC_Point, EC_Point);

  // EC_Double - double a point on the curve
  EC_Point EC_Double(EC_Point);

  // EC_Basic_Multiply - brute-force point multiplication
  EC_Point EC_Basic_Multiply(EC_Point, BigInteger);

  // EC_Double_and_Add - exponential point multiplication
  EC_Point EC_Double_and_Add(EC_Point, BigInteger);

  // EC_WindowNAF - window-NAF method point multiplication - improved
  // performance
  EC_Point EC_WindowNAF(EC_Point, BigInteger, int32_t);

  // Helper functions
  void EC_Find_NAF(int32_t, std::vector<int32_t> &, std::vector<int32_t> &);
  EC_Point EC_Negative(EC_Point &);
};
