#ifndef BIGINTEGER_H // Include guard
#define BIGINTEGER_H

#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<sys/time.h>
#include<climits>
#include<cstdlib>

const int BI_BASE = 1000000;
const int BI_BASE_DIGITS = 6;

// forward declarations:
class BigInteger;
bool BI_Miller_Rabin( BigInteger p, int it = 32 );

class BigInteger {

    // Stores arbitrarily large integer absolute value
    std::vector<int> a;
    int sign;

public:

// Helper functions
    void read( const std::string & );
    void trim();
    friend std::vector<long long> karatsuba( const std::vector<long long>&, const std::vector<long long>& );
    friend std::vector<int> BI_Helper_DecBin( const BigInteger& );
    friend BigInteger BI_Helper_BinDec( const std::vector<int>& );
    friend std::vector<int> BI_Helper_Subtract( const std::vector<int>&, const std::vector<int>& );
    friend bool BI_Helper_Less( std::vector<int>, std::vector<int> );
    friend int BI_Helper_CompareAbs( const BigInteger&, const BigInteger& );
    friend std::pair<BigInteger, BigInteger> BI_Helper_Divide( const BigInteger&, const BigInteger& );

// Constructors
    BigInteger();
    BigInteger( int );
    BigInteger( long long );
    BigInteger( const BigInteger& );
    BigInteger( const std::string& );

// operator=
    BigInteger& operator=( int );
    BigInteger& operator=( long long );
    BigInteger& operator=( const BigInteger& );

// istream and ostream operators
    friend std::istream& operator>>( std::istream&, BigInteger& );
    friend std::ostream& operator<<( std::ostream&, const BigInteger& );

// Relational operators
    bool operator<( const BigInteger& ) const;
    bool operator>( const BigInteger& ) const;
    bool operator<=( const BigInteger& ) const;
    bool operator>=( const BigInteger& ) const;
    bool operator==( const BigInteger& ) const;
    bool operator!=( const BigInteger& ) const;

// long long relational operators
    bool operator<( const long long ) const;
    bool operator>( const long long ) const;
    bool operator<=( const long long ) const;
    bool operator>=( const long long ) const;
    bool operator==( const long long ) const;
    bool operator!=( const long long ) const;

// Compound arithmetic operators
    BigInteger& operator+=( const BigInteger& );
    BigInteger& operator-=( const BigInteger& );
    BigInteger& operator*=( const BigInteger& );
    BigInteger& operator/=( const BigInteger& );
    BigInteger& operator%=( const BigInteger& );

// Arithmetic operators (non-members)
    friend const BigInteger operator+( BigInteger, const BigInteger& );
    friend const BigInteger operator-( BigInteger, const BigInteger& );
    friend const BigInteger operator*( BigInteger, const BigInteger& );
    friend const BigInteger operator/( BigInteger, const BigInteger& );
    friend const BigInteger operator%( BigInteger, const BigInteger& );

// long long compound arithmetic operators
    BigInteger& operator+=( long long );
    BigInteger& operator-=( long long );
    BigInteger& operator*=( long long );
    BigInteger& operator/=( long long );
    BigInteger& operator%=( long long );

// long long arithmetic operators (non-members)
    friend const BigInteger operator+( BigInteger, long long );
    friend const BigInteger operator-( BigInteger, long long );
    friend const BigInteger operator*( BigInteger, long long );
    friend const BigInteger operator/( BigInteger, long long );
    friend const BigInteger operator%( BigInteger, long long );

// pre- and post-increment and decrement operators
    BigInteger& operator++();
    BigInteger operator++( int );
    BigInteger& operator--();
    BigInteger operator--( int );

// Miscellaneous
    long long to_llong() const;
    bool is_even() const;
    bool is_zero() const;
    bool is_negative() const;
    friend BigInteger abs( BigInteger& );
    friend BigInteger gcd( BigInteger, BigInteger );
    friend BigInteger lcm( BigInteger, BigInteger );
    friend BigInteger BI_FastExp( BigInteger, BigInteger, BigInteger );
    friend BigInteger BI_FastExp( BigInteger, long long, BigInteger );
    friend int BI_FastExp( BigInteger, long long, int );
    friend BigInteger BI_ModInv( BigInteger, BigInteger );
    friend BigInteger BI_Residue( BigInteger, BigInteger );
    friend bool BI_Miller_Rabin( BigInteger, int );
};

#endif
