#include "BigInteger.hpp"
#include<iomanip>
#include<algorithm>
#include<cstdio>
#include<sys/time.h>
#include<climits>

using namespace std;

// Helper functions

// read - read from string
void BigInteger::read( const string& v ) {
    a.clear();
    sign = 1;
    int i=0, j, k, l, digit, flag = 0;
    while( i < (int)(v.size()) && ( v[i] == '+' || v[i] == '-' || v[i] == ' ' ) ) {
        if( v[i] == '-' ) sign = -sign;
        ++i;
    }
    j = i;
    while( j < (int)(v.size()) && v[j] >= '0' && v[j] <= '9' ) ++j;
    --j;
    if( j < i ) {
        a.push_back(0);
        sign = 1;
        return;
    }
    for( ; j >= i; j -= BI_BASE_DIGITS ) {
        k = 0;
        for( l = max( i, j+1-BI_BASE_DIGITS ); l <= j; ++l ) {
            k = k*10 + v[l] - '0';
        }
        a.push_back(k);
    }
    trim();
}

// trim - trim leading zeros
void BigInteger::trim() {
    while ( !a.empty() && !a.back() ) a.pop_back();
    if( a.empty() ) sign = 1;
}

// karatsuba - fast multiplication
vector<long long> karatsuba( const vector<long long>& A, const vector<long long>& B ) {
    int i, j, n = A.size(), m;
    vector<long long> res( n+n, 0 );
    if( n <= 512 ) { // Switch to basic (close to optimal for ~32 decimal digits)
        for( i = 0; i < n; ++i ) {
            for( j = 0; j < n; ++j ) res[i+j] += A[i] * B[j];
        }
        return res;
    }

    m = n >> 1;
    vector<long long> a1( A.begin(), A.begin() + m );
    vector<long long> a2( A.begin() + m, A.end() );
    vector<long long> b1( B.begin(), B.begin() + m );
    vector<long long> b2( B.begin() + m, B.end() );

    vector<long long> x = karatsuba( a1, b1 );
    vector<long long> y = karatsuba( a2, b2 );
    for( i = 0; i < m; ++i ) {
        a2[i] += a1[i];
        b2[i] += b1[i];
    }
    vector<long long> r = karatsuba( a2, b2 );

    for( i = 0; i < x.size(); ++i ) r[i] -= x[i];
    for( i = 0; i < y.size(); ++i ) r[i] -= y[i];
    for( i = 0; i < r.size(); ++i ) res[i+m] += r[i];
    for( i = 0; i < x.size(); ++i ) res[i] += x[i];
    for( i = 0; i < y.size(); ++i ) res[i] += y[i];
    return res;
}

// BI_Helper_DecBin - convert to binary
vector<int> BI_Helper_DecBin( const BigInteger& bi ) {
    int i, j, x;
    vector<int> dbi, temp;
    for( i = bi.a.size()-1; i >= 0; --i ) {
        x = bi.a[i];
        vector<int> p;
        while(x) {
            p.push_back(x%10);
            x /= 10;
        }
        while( p.size() < BI_BASE_DIGITS ) p.push_back(0);
        for( j = p.size()-1; j >= 0; --j ) dbi.push_back(p[j]);
    }
    j = 0; while( j < dbi.size() && dbi[j] == 0 ) ++j;
    dbi.erase( dbi.begin(), dbi.begin() + j );

    vector<int> bin;
    while(dbi.size()) {
        bin.push_back( dbi[dbi.size()-1] & 1 );
        temp.clear();
        int adt = 0, dig, flag = 0;
        for( i = 0; i < dbi.size(); ++i ) {
            dig = (dbi[i]/2) + adt;
            if( dig ) flag = 1;
            if( flag ) temp.push_back(dig);
            adt = (dbi[i]&1) ? 5 : 0;
        }
        dbi.swap(temp);
    }
    reverse( bin.begin(), bin.end() );
    return bin;
}

// BI_Helper_BinDec - convert from binary
BigInteger BI_Helper_BinDec( const vector<int>& vi ) {
    BigInteger b(1), res(0);
    for( int i = vi.size()-1; i >= 0; --i ) {
        if( vi[i] == 1 ) {
            res += b;
        }
        b += b;
    }
    return res;
}

// BI_Helper_Subtract - subtract two numbers in binary
vector<int> BI_Helper_Subtract( const vector<int>& a, const vector<int>& b ) {
    // find a-b in binary representation, assumes a > b
    int i, j;
    for( i = 0; i < b.size(); ++i ) if( b[i] ) break;
    if( i == b.size() ) { // b is zero
        return a;
    }
    vector<int> bc; // two's complement of b
    for( int i = b.size(); i < a.size(); ++i ) bc.push_back(1);
    for( i = 0; i < b.size(); ++i ) bc.push_back(1 - b[i]);
    int d, c = 1; // Initial carry of 1 converts 1's complement to 2's complement

    vector<int> res(a.size(), 0);
    for( int i = a.size()-1; i >= 0; --i ) {
        d = c + a[i] + bc[i];
        res[i] = (d & 1); // d % 2
        c = (d >> 1); // d / 2
    }

    j = 0; while( j < res.size() && res[j] == 0 ) ++j;
    res.erase( res.begin(), res.begin() + j );
    return res;
}

// BI_Helper_Less - similar to std::Less for BigInteger objects
bool BI_Helper_Less( vector<int> a, vector<int> b ) {
    int j = 0;
    while( j < a.size() && a[j] == 0 ) ++j;
    a.erase( a.begin(), a.begin() + j );
    j = 0;
    while( j < b.size() && b[j] == 0 ) ++j;
    b.erase( b.begin(), b.begin() + j );

    if( a.size() < b.size() ) return true;
    if( b.size() < a.size() ) return false;
    for( j = 0; j < a.size(); ++j ) {
        if( a[j] == 1 && b[j] == 0 ) return false;
        if( a[j] == 0 && b[j] == 1 ) return true;
    }
    return false;
}

// BI_Helper_CompareAbs - compare absolute values - 0 equal, 1 greater, 2 less
int BI_Helper_CompareAbs( const BigInteger& a1, const BigInteger& b1 ) {
    // this loop runs in approx same time always as desired
    int sflag = 0;
    if( a1.a.size() > b1.a.size() ) sflag = 1;
    if( a1.a.size() < b1.a.size() ) sflag = 2;
    else if( a1.a.size() == b1.a.size() ) {
        for( int i = a1.a.size() - 1; i >= 0; --i ) {
            if( a1.a[i] != b1.a[i] && sflag == 0 ) {
                if( a1.a[i] > b1.a[i] ) sflag = 1;
                else sflag = 2;
            }
        }
    }
    return sflag;
}

// BI_Helper_Divide - find quotient and remainder
pair<BigInteger, BigInteger> BI_Helper_Divide( const BigInteger& a1, const BigInteger& b1 ) {
    int i, j, x, y;
    int sflag = BI_Helper_CompareAbs(a1,b1);
    switch(sflag) {
        case 0 : // equal - also handles 0 / 0
            {
                BigInteger xx(1), yy(0);
                if( a1.sign != b1.sign ) xx.sign = -1;
                if( a1.a.size() == 0 ) return pair<BigInteger,BigInteger>(yy,yy);
                return pair<BigInteger,BigInteger>(xx,yy);
            }
            break;
        case 1 : // abs(a1) > abs(b1) - also handles a1 / 0
            {
                // Convert both absolute values to base 2 vectors
                    vector<int> a = BI_Helper_DecBin(a1);
                    vector<int> b = BI_Helper_DecBin(b1);

                // Divide a by b;
                    if( b.empty() ) {
                        BigInteger xx(0), yy(0);
                        return make_pair(xx,yy);
                    }
                    vector<int> r, q;
                    i = 0;
                    while( 1 ) {
                        int t = 0;
                        while( t < r.size() && r[t] == 0 ) ++t;
                        r.erase( r.begin(), r.begin() + t );

                        if( BI_Helper_Less(r,b) ) {
                            q.push_back(0);
                        }
                        else {
                            r = BI_Helper_Subtract(r,b);
                            q.push_back(1);
                        }
                        if( i < a.size() ) r.push_back(a[i++]);
                        else break;
                    }
                BigInteger xx = BI_Helper_BinDec(q);
                BigInteger yy = BI_Helper_BinDec(r);
                xx.sign = a1.sign * b1.sign;
                yy.sign = a1.sign;
                xx.trim(); yy.trim();
                return pair<BigInteger,BigInteger>(xx,yy);
            }
            break;
        case 2 : // abs(a1) < abs(b1)
            {
                BigInteger xx(0), yy(a1);
                xx.trim(); yy.trim();
                return pair<BigInteger,BigInteger>(xx,yy);
            }
            break;
    }
    BigInteger xx(0), yy(0);
    return pair<BigInteger,BigInteger>(xx,yy);
}

// Constructors

BigInteger::BigInteger() {
    a.clear();
    sign = 1;
}
BigInteger::BigInteger( int v ) {
    *this = v;
}
BigInteger::BigInteger( long long v ) {
    *this = v;
}
BigInteger::BigInteger( const BigInteger& v ) {
    *this = v;
}
BigInteger::BigInteger( const string& v ) {
    read(v);
}

// operator=

BigInteger& BigInteger::operator=( int v ) {
    a.clear();
    if( v < 0 ) {
        v = -v;
        sign = -1;
    }
    else sign = 1;
    for( ; v > 0; v /= BI_BASE ) {
        a.push_back( v%BI_BASE );
    }
    return *this;
}
BigInteger& BigInteger::operator=( long long v ) {
    a.clear();
    if( v < 0 ) {
        v = -v;
        sign = -1;
    }
    else sign = 1;
    for( ; v > 0; v /= BI_BASE ) {
        a.push_back( v%BI_BASE );
    }
    return *this;
}
BigInteger& BigInteger::operator=( const BigInteger& v ) {
    if( this != &v ) {
        a = v.a;
        sign = v.sign;
    }
    trim();
    return *this;
}

// istream and ostream operators

istream& operator>>( istream& stream, BigInteger& v ) {
    string s;
    stream >> s;
    v.read(s);
    return stream;
}
ostream& operator<<( ostream& stream, const BigInteger& v ) {
    if( v.sign == -1 ) stream << '-';
    stream << ( v.a.empty() ? 0 : v.a.back() );
    for( int i = v.a.size() - 2; i >= 0; --i ) {
        stream << setw(BI_BASE_DIGITS) << setfill('0') << v.a[i];
    }
    return stream;
}

// Relational operators

bool BigInteger::operator<( const BigInteger& y ) const {
    if( sign != y.sign ) return sign < y.sign;
    if( a.size() != y.a.size() ) return a.size()*sign < y.a.size()*y.sign;
    for( int i = a.size()-1; i >= 0; --i ) {
        if( a[i] != y.a[i] ) return a[i]*sign < y.a[i]*sign;
    }
    return false;
}
bool BigInteger::operator>( const BigInteger& y ) const {
    return y < *this;
}
bool BigInteger::operator<=( const BigInteger& y ) const {
    return !( *this > y );
}
bool BigInteger::operator>=( const BigInteger& y ) const {
    return !( *this < y );
}
bool BigInteger::operator==( const BigInteger& y ) const {
    if( sign != y.sign ) return false;
    if( a.size() != y.a.size() ) return false;
    for( int i = 0; i < a.size(); ++i ) {
        if( a[i] != y.a[i] ) return false;
    }
    return true;
}
bool BigInteger::operator!=( const BigInteger& y ) const {
    return !( *this == y );
}

// long long relational operators

bool BigInteger::operator<( long long rhs ) const {
    BigInteger y(rhs);
    if( sign != y.sign ) return sign < y.sign;
    if( a.size() != y.a.size() ) return a.size()*sign < y.a.size()*y.sign;
    for( int i = a.size()-1; i >= 0; --i ) {
        if( a[i] != y.a[i] ) return a[i]*sign < y.a[i]*sign;
    }
    return false;
}
bool BigInteger::operator>( long long rhs ) const {
    BigInteger y(rhs);
    return y < *this;
}
bool BigInteger::operator<=( long long rhs ) const {
    BigInteger y(rhs);
    return !( *this > y );
}
bool BigInteger::operator>=( long long rhs ) const {
    BigInteger y(rhs);
    return !( *this < y );
}
bool BigInteger::operator==( long long rhs ) const {
    BigInteger y(rhs);
    if( sign != y.sign ) return false;
    if( a.size() != y.a.size() ) return false;
    for( int i = 0; i < a.size(); ++i ) {
        if( a[i] != y.a[i] ) return false;
    }
    return true;
}
bool BigInteger::operator!=( long long rhs ) const {
    BigInteger y(rhs);
    return !( *this == y );
}

// Compound arithmetic operators

BigInteger& BigInteger::operator+=( const BigInteger& rhs ) {
    BigInteger y(rhs);
    if( sign == y.sign ) {
        // add absolute values
        for( int i = 0, c = 0; ( i < (int)(max(a.size(), y.a.size())) ) || c>0; ++i ) {
            if( i == (int)(a.size()) ) a.push_back(0);
            a[i] += c + ( i<(int)(y.a.size()) ? y.a[i] : 0 );
            if( a[i] >= BI_BASE ) {
                c = 1;
                a[i] -= BI_BASE;
            }
            else c = 0;
        }
    }
    else {
        int sflag = BI_Helper_CompareAbs(*this,y);
        if( sflag == 1 ) {
            // abs(this) > abs(y), this = this - y and keep sign(this)
            for( int i = 0, c = 0; i < y.a.size() || c > 0; ++i ) {
                a[i] -= c + ( i < y.a.size() ? y.a[i] : 0 );
                if( a[i] < 0 ) {
                    c = 1;
                    a[i] += BI_BASE;
                }
                else c = 0;
            }
        }
        else {
            // abs(this) <= abs(y), this = y - this and keep sign(y)
            BigInteger t(y);
            for( int i = 0, c = 0; i < a.size() || c > 0; ++i ) {
                t.a[i] -= c + ( i < a.size() ? a[i] : 0 );
                if( t.a[i] < 0 ) {
                    c = 1;
                    t.a[i] += BI_BASE;
                }
                else c = 0;
            }
            *this = t;
        }
    }
    trim();
    return *this;
}
BigInteger& BigInteger::operator-=( const BigInteger& rhs ) {
    BigInteger y(rhs);
    if( sign != y.sign ) {
        // add absolute values
        for( int i = 0, c = 0; ( i < (int)(max(a.size(), y.a.size())) ) || c>0; ++i ) {
            if( i == (int)(a.size()) ) a.push_back(0);
            a[i] += c + ( i<(int)(y.a.size()) ? y.a[i] : 0 );
            if( a[i] >= BI_BASE ) {
                c = 1;
                a[i] -= BI_BASE;
            }
            else c = 0;
        }
    }
    else {
        // compare abs(this) and abs(y) - note this loop runs in approx same time always as desired
        int sflag = BI_Helper_CompareAbs(*this,y);
        if( sflag == 1 ) {
            // abs(this) > abs(y), this = this - y and keep sign(this)
            for( int i = 0, c = 0; i < y.a.size() || c > 0; ++i ) {
                a[i] -= c + ( i < y.a.size() ? y.a[i] : 0 );
                if( a[i] < 0 ) {
                    c = 1;
                    a[i] += BI_BASE;
                }
                else c = 0;
            }
        }
        else {
            // abs(this) <= abs(y), this = y - this and flip sign(y)
            BigInteger t(y);
            for( int i = 0, c = 0; i < a.size() || c > 0; ++i ) {
                t.a[i] -= c + ( i < a.size() ? a[i] : 0 );
                if( t.a[i] < 0 ) {
                    c = 1;
                    t.a[i] += BI_BASE;
                }
                else c = 0;
            }
            *this = t;
            sign = -sign;
        }
    }
    trim();
    return *this;
}
BigInteger& BigInteger::operator*=( const BigInteger& rhs ) {
    BigInteger y(rhs);
    vector<long long> v1( a.begin(), a.end() );
    vector<long long> v2( y.a.begin(), y.a.end() );

    // Pad v1 and v2 with zeros until their length is equal and a multiple of 2
    int s1 = v1.size(), s2 = v2.size();
    while( s1 < s2 ) v1.push_back(0), ++s1;
    while( s2 < s1 ) v2.push_back(0), ++s2;
    while( v1.size() & (v1.size() - 1) ) v1.push_back(0), v2.push_back(0);

    vector<long long> c = karatsuba(v1, v2);
    a.clear();
    sign *= y.sign;
    for( int i = 0, carry = 0; i < c.size(); ++i ) {
        long long cur = c[i] + carry;
        a.push_back( (int)( cur%BI_BASE ) );
        carry = (int)( cur/BI_BASE );
    }
    trim();
    return *this;
}
BigInteger& BigInteger::operator/=( const BigInteger& rhs ) {
    BigInteger y(rhs);
    pair<BigInteger, BigInteger> result = BI_Helper_Divide(*this, y);
    *this = result.first;
    return *this;
}
BigInteger& BigInteger::operator%=( const BigInteger& rhs ) {
    BigInteger y(rhs);
    pair<BigInteger, BigInteger> result = BI_Helper_Divide(*this, y);
    *this = result.second;
    return *this;
}

// Arithmetic operators (non-members)

const BigInteger operator+( BigInteger x, const BigInteger& y ) {
    x += y;
    return x;
}
const BigInteger operator-( BigInteger x, const BigInteger& y ) {
    x -= y;
    return x;
}
const BigInteger operator*( BigInteger x, const BigInteger& y ) {
    x *= y;
    return x;
}
const BigInteger operator/( BigInteger x, const BigInteger& y ) {
    x /= y;
    return x;
}
const BigInteger operator%( BigInteger x, const BigInteger& y ) {
    x %= y;
    return x;
}

// long long compound arithmetic operators

BigInteger& BigInteger::operator+=( long long y ) {
    *this += BigInteger(y);
    return *this;
}
BigInteger& BigInteger::operator-=( long long y ) {
    *this -= BigInteger(y);
    return *this;
}
BigInteger& BigInteger::operator*=( long long y ) {
    *this *= BigInteger(y);
    return *this;
}
BigInteger& BigInteger::operator/=( long long y ) {
    *this /= BigInteger(y);
    return *this;
}
BigInteger& BigInteger::operator%=( long long y ) {
    *this %= BigInteger(y);
    return *this;
}

// long long arithmetic operators (non-members)

const BigInteger operator+( BigInteger x, long long y ) {
    x += y;
    return x;
}
const BigInteger operator-( BigInteger x, long long y ) {
    x -= y;
    return x;
}
const BigInteger operator*( BigInteger x, long long y ) {
    x *= y;
    return x;
}
const BigInteger operator/( BigInteger x, long long y ) {
    x /= y;
    return x;
}
const BigInteger operator%( BigInteger x, long long y ) {
    x %= y;
    return x;
}

// pre- and post-increment and decrement operators

BigInteger& BigInteger::operator++() {
    trim();
    int i = 0;
    if( sign == 1 ) { // handles 0
        while( i < a.size() && (++a[i]) == BI_BASE ) a[i++] = 0;
        if( i == a.size() ) a.push_back(1);
    }
    else {
        while( i < a.size() && (--a[i]) == -1 ) a[i++] = BI_BASE-1;
    }
    trim();
    return *this;
}
BigInteger BigInteger::operator++(int) {
    BigInteger temp(*this);
    operator++();
    return temp;
}
BigInteger& BigInteger::operator--() {
    trim();
    int i = 0;
    if( a.size() == 0 ) {
        a.push_back(1);
        sign = -1;
    }
    else if( sign == 1 ) {
        while( i < a.size() && (--a[i]) == -1 ) a[i++] = BI_BASE-1;
    }
    else {
        while( i < a.size() && (++a[i]) == BI_BASE ) a[i++] = 0;
        if( i == a.size() ) a.push_back(1);
    }
    trim();
    return *this;
}
BigInteger BigInteger::operator--(int) {
    BigInteger temp(*this);
    operator--();
    return temp;
}

// Miscellaneous

// to_llong - returns a long long int with the value of (*this). Value undefined if outside LLONG limits
long long BigInteger::to_llong() const {
    long long res = 0;
    for( int i = a.size() - 1; i >= 0; i-- ) res = res * BI_BASE + a[i];
    return res * sign;
}

// is_even - returns true if even
bool BigInteger::is_even() const {
    if( a.size() == 0 ) return true;
    if( a[0] & 1 ) return false;
    return true;
}

// is_zero - returns true if zero
bool BigInteger::is_zero() const {
    if( a.size() == 0 ) return true;
    return false;
}

// is_negative - returns true if negative
bool BigInteger::is_negative() const {
    if( sign == -1 ) return true;
    return false;
}

// abs - returns absolute value of x
BigInteger abs( BigInteger& x ) {
    BigInteger res(x);
    res.sign = 1;
    return res;
}

// gcd - returns greatest common positive divisor
BigInteger gcd( BigInteger m, BigInteger n ) {
    if( m.sign == -1 ) m = abs(m);
    if( n.sign == -1 ) n = abs(n);
    if( n == 0 ) return m;
    return gcd( n, m % n );
}

// lcm - returns least common positive multiple
BigInteger lcm( BigInteger m, BigInteger n ) {
    if( m.sign == -1 ) m = abs(m);
    if( n.sign == -1 ) n = abs(n);
    return ( m * n ) / gcd( m, n );
}

// BI_FastExp - modular exponentiation by repeated squaring
BigInteger BI_FastExp( BigInteger base, BigInteger exp, BigInteger m ) {
    m.trim(); if( m.sign == -1 || m.a.size() == 0 || exp < 0 ) return *(new BigInteger());
    base.trim(); exp.trim();
    if( exp.is_zero() ) return *(new BigInteger(1));
    vector<int> e = BI_Helper_DecBin(exp);
    BigInteger res = 1;
    for( int i = e.size()-1; i >= 0; --i ) {
        if( e[i] ) res = (res * base) % m;
        base = ( base * base ) % m;
    }
    res.trim();
    return res;
}
BigInteger BI_FastExp( BigInteger base, long long exp, BigInteger m ) {
    m.trim(); if( m.sign == -1 || m.a.size() == 0 || exp < 0 ) return *(new BigInteger());
    base.trim();
    BigInteger res = 1;
    while( exp ) {
        if( exp & 1 ) res = (res * base) % m;
        base = ( base * base ) % m;
        exp >>= 1;
    }
    res.trim();
    return res;
}

// BI_FastExp - modular exponentiation by repeated squaring - small modulus
int BI_FastExp( BigInteger base, long long exp, int m ) {
    if( m <= 0 || exp < 0 ) return 0;
    base.trim(); base %= m;
    long long my_base = base.to_llong();
    long long res = 1;
    while( exp ) {
        if( exp & 1 ) res = (res * my_base) % m;
        my_base = ( my_base * my_base ) % m;
        exp >>= 1;
    }
    return (int)(res%m);
}

// BI_ModInv - modular inverse
BigInteger BI_ModInv( BigInteger a, BigInteger m ) {
    m.trim();
    if( m.sign == -1 || m.a.size() == 0 ) return *(new BigInteger());
    BigInteger m0(m), x0(0), x1(1);
    BigInteger t, q;
	if( m == 1 ) return 1;
	while( a > 1 ) {
		q = a / m;
		t = m;
		m = a % m;
		a = t;
		t = x0;
		x0 = x1 - (q*x0);
		x1 = t;
	}
	if( x1 < 0 ) x1 += m0;
	x1.trim();
	return x1;
}

// BI_Residue - Returns least nonnegative residue of a modulo m
BigInteger BI_Residue( BigInteger a, BigInteger m ) {
    m.trim();
    if( m.sign == -1 || m.a.size() == 0 ) return *(new BigInteger());
    BigInteger res = a % m;
    if( res.sign == -1 ) return res + m;
    return res;
}

// BI_Miller_Rabin - Returns true if p is prime, false with high probability if p is composite
bool BI_Miller_Rabin( BigInteger p, int it ) {
    if( p < 2 ) return false;
    if( p <= 10 ) {
        if( p == 2 || p == 3 || p == 5 || p == 7 ) return true;
        else return false;
    }
    if( p.is_even() ) return false;

    srand(time(NULL));

    int k = 0;
    BigInteger d(p);
    --d;
    while( d.is_even() ) {
        d /= 2;
        ++k;
    }

    // Security level and performance guards
    if( it < 32 ) it = 32;
    if( it > 128 ) it = 128;
    long long rmax = RAND_MAX, ran;
    int rm;
    if( p >= rmax ) rm = rmax - 2;
    else {
        BigInteger tp(p);
        --tp; --tp;
        rm = (int)(tp.to_llong());
    }

    for( long long i = 0; i < it; ++i ) {
        ran = (rand()%rm) + 2;
        BigInteger a = BigInteger( ran );
        BigInteger val = BI_FastExp( a, d, p );
        if( val % p == 1 ) {
            continue; // Not a witness
        }
        BigInteger chk = p-1;
        int flag = 0;
        for( int g = 0; g < k; ++g ) {
            if( val == chk ) {
                flag = 1; break; // Not a witness
            }
            val = BI_Residue( val * val, p );
        }
        if( flag ) continue;
        // Is a witness for compositeness
        return false; // Return composite
    }
    return true; // Return probable prime
}


