#include<iomanip>
#include<algorithm>
#include<cstdio>
#include<sys/time.h>
#include<climits>
#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include "BigInteger.hpp"
#include "ECPO.hpp"

using namespace std;

// EC_Point functions

// Constructors

EC_Point::EC_Point() {
    zero = false;
    x = BigInteger(0);
    y = BigInteger(0);
}

EC_Point::EC_Point( int z ) {
    if( z == 0 ) zero = true;
    else zero = false;
    x = BigInteger(0);
    y = BigInteger(0);
}

EC_Point::EC_Point( const BigInteger& xx, const BigInteger& yy ) {
    x = xx;
    y = yy;
    zero = false;
}

EC_Point::EC_Point( const EC_Point& p1 ) {
    zero = p1.zero;
    x = p1.x;
    y = p1.y;
}

EC_Point::EC_Point( pair<BigInteger,BigInteger>& P ) {
    x = P.first;
    y = P.second;
    zero = false;
}

// ostream operator

ostream& operator<<( ostream& stream, const EC_Point& pt ) {
    if( pt.zero ) {
        stream << "O";
        return stream;
    }
    stream << pt.x;
    stream << ", ";
    stream << pt.y;
    return stream;
}

// EC_Curve functions

// Constructors

EC_Curve::EC_Curve() {
    p = 0;
    A = 0;
    B = 0;
    valid = false;
    warning = false;
}
EC_Curve::EC_Curve( const BigInteger& prime, const BigInteger& Aval, const BigInteger& Bval ) {
    // change to least residues here
    p = prime;
    A = Aval;
    B = Bval;
    valid = false;
    warning = false;
}

// EC_Validate - validate the curve for point operations (p prime, discriminant not 0)

void EC_Curve::EC_Validate() {
    if( (A*A*A*4) + (B*B*27) != 0 ) {
        if( BI_Miller_Rabin( p ) ) valid = true;
        else {
            cerr << "WARNING : Curve parameter not prime" << endl;
            valid = false;
        }
    }
    else {
        cerr << "WARNING : Curve discriminant zero" << endl;
        valid = false;
    }
    warning = false;
}

// EC_Validate_Point - returns true if the parameter point is on the curve

bool EC_Curve::EC_Validate_Point( EC_Point pt ) {
    if( !valid && !warning ) {
        std::cerr << "WARNING : Curve invalid / not validated" << std::endl;
        warning = true;
    }
    if( pt.zero ) return true;
    if( pt.x >= p || pt.x < 0 || pt.y >= p || pt.y < 0 ) return false;
    BigInteger left = BI_FastExp( pt.y, 2, p );
    BigInteger right = BI_FastExp( pt.x, 3, p ) + ( (pt.x*A) % p ) + B;
    left = BI_Residue( left, p );
    right = BI_Residue( right, p );
    return left == right;
}

// EC_Add - add two points on the curve

EC_Point EC_Curve::EC_Add( EC_Point p1, EC_Point p2 ) {
    if( !valid && !warning ) {
        std::cerr << "WARNING : Curve invalid / not validated" << std::endl;
        warning = true;
    }
    if( !EC_Validate_Point(p1) ) {
        std::cerr << "Error(2) : point not on curve - " << std::endl << p1.x << std::endl << p1.y << std::endl;
        return *(new EC_Point(0));
    }
    if( !EC_Validate_Point(p2) ) {
        std::cerr << "Error(2) : point not on curve - " << std::endl << p2.x << std::endl << p2.y << std::endl;
        return *(new EC_Point(0));
    }

    if( p1.zero ) return p2;
    if( p2.zero ) return p1;
    BigInteger temp = BI_Residue( p1.y + p2.y, p );
    if( p1.x == p2.x && temp == 0 ) return *(new EC_Point(0));

    BigInteger lam_n, lam_d, lam;
    if( p1.x != p2.x ) {
        lam_n = BI_Residue( p2.y - p1.y, p );
        lam_d = BI_Residue( p2.x - p1.x, p );
        lam_d = BI_ModInv( lam_d, p );
        lam = BI_Residue( lam_n * lam_d, p );
    }
    else if( p1.x == p2.x && p1.y == p2.y ){
        BigInteger t = BI_Residue( p1.x * p1.x * 3, p );
        lam_n = BI_Residue( t + A, p );
        lam_d = BI_Residue( p1.y * 2, p );
        lam_d = BI_ModInv( lam_d, p );
        lam = BI_Residue( lam_n * lam_d, p );
    }
    else {
        std::cerr << "Error(0) : something went wrong" << std::endl;
        return *(new EC_Point(0));
    }

    EC_Point res;
    res.x = BI_Residue( lam*lam, p ) - p1.x - p2.x;
    res.x = BI_Residue( res.x, p );
    res.y = lam * BI_Residue( p1.x - res.x, p ) - p1.y;
    res.y = BI_Residue( res.y, p );

    if( EC_Validate_Point(res) ) return res;
    else {
        std::cerr << "Error(1) : something went wrong" << std::endl;
        return *(new EC_Point(0));
    }
}

// EC_Double - double a point on the curve

EC_Point EC_Curve::EC_Double( EC_Point p1 ) {
    return EC_Add(p1, p1);
}

// EC_Basic_Multiply - brute-force point multiplication

EC_Point EC_Curve::EC_Basic_Multiply( EC_Point p1, BigInteger n ) {
    if( n <= 0 ) return *(new EC_Point(0));
    EC_Point res(0);
    while( n.is_zero() == false ) {
        res = EC_Add( res, p1 );
        --n;
    }
    return res;
}

// EC_Double_and_Add - exponential point multiplication

EC_Point EC_Curve::EC_Double_and_Add( EC_Point p1, BigInteger n ) {
    if( n <= 0 ) return *(new EC_Point(0));
    vector<int> e = BI_Helper_DecBin(n);
    EC_Point res(0);
    for( int i = e.size()-1; i >= 0; --i ) {
        if( e[i] ) res = EC_Add( res, p1 );
        p1 = EC_Double( p1 );
    }
    return res;
}

// EC_WindowNAF - window-NAF method point multiplication - improved performance
// Algorithm 3.36 from Guide to Elliptic Curve Cryptography

EC_Point EC_Curve::EC_WindowNAF( EC_Point p1, BigInteger n, int w ) {
    if( n <= 0 ) return *(new EC_Point(0));
    if( n <= (1<<w) ) return EC_Double_and_Add(p1,n);
    vector<int> k = BI_Helper_DecBin(n);
    vector<int> naf;

    // Correctness and performance guards

        if( w < 2 ) w = 2;
        if( w > 16 ) w = 16;

    // Compute non-adjacent form of the binary representation of n in an integer vector

        EC_Find_NAF( w, k, naf);

    // Compute 2^i P1 for i in { 1, 3, 5, ..., lim }, and their negatives.

        int lim = (1<<(w-1)) - 1, i;

        EC_Point dbl = EC_Double(p1);
        vector< EC_Point > vec;
        vector< EC_Point > neg;
        vec.reserve(lim+1);
        neg.reserve(lim+1);

        vec.push_back(*(new EC_Point(0)));
        vec.push_back(p1);
        for( i = 3; i <= lim; i += 2 ) {
            vec.push_back(*(new EC_Point(0)));
            p1 = EC_Add(p1,dbl);
            vec.push_back(p1);
        }
        for( i = 0; i <= lim; ++i ) neg.push_back(EC_Negative(vec[i]));

    // Compute the point multiple

        EC_Point res(0);
        int nsz = naf.size();
        for( i = 0; i < nsz; ++i ) {
            res = EC_Double(res);
            if( naf[i] > 0 ) {
                res = EC_Add( res, vec[naf[i]] );
            }
            else if( naf[i] < 0 ) {
                res = EC_Add( res, neg[-naf[i]] );
            }
        }
        return res;

}

// Helper functions

// EC_Find_NAF - constructs the width-w non-adjacent form vector of k in naf
// Algorithm 3.35 from Guide to Elliptic Curve Cryptography

void EC_Curve::EC_Find_NAF( int w, vector<int>& k, vector<int>& naf ) {
    naf.clear();
    int sz = k.size();
    naf.resize(sz+1, 0);
    reverse( k.begin(), k.end() );
    k.resize(sz+w+1,0);
    int msk, i, j, ct;

    int num = (1<<w) - 1, cmp = (1<<(w-1)), x;
    for( i = 0; i < sz; ++i ) {
        if( k[i] ) {
            msk = 0; ct = 0;
            for( j = i; ct < w; ++j ) {
                msk |= ((k[j])<<(ct++));
            }
            if( msk >= cmp ) { // negative
                msk = (msk^num) + 1;
                naf[i] = -msk;
                for( x = 1; x < w; ++x ) k[i+x] = 0;
                for( x = j; k[x] == 1; ++x ) {
                    k[x] = 0;
                }
                k[x] = 1;
            }
            else {
                naf[i] = msk;
                for( x = 1; x < w; ++x ) k[i+x] = 0;
            }
        }
        else {
            naf[i] = 0;
        }
    }
    while( naf.back() == 0 ) naf.pop_back();
    reverse( naf.begin(), naf.end() );
}

// EC_Negative - return negative of a point on the curve

EC_Point EC_Curve::EC_Negative( EC_Point& p1 ) {
    EC_Point p2(p1);
    p2.y = p - p2.y;
    return p2;
}




