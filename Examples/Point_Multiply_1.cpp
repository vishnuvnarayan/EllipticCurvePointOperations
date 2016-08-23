#include<iostream>
#include<sys/time.h>
#include<cstdlib>
#include "BigInteger.hpp"
#include "ECPO.hpp"
using namespace std;

int main() {
    BigInteger p("13"), A("3"), B("8");

    EC_Curve C(p,A,B);
    C.EC_Validate();

    EC_Point P1( BigInteger("1"), BigInteger("5") );
    EC_Point P2( BigInteger("2"), BigInteger("3") );

    BigInteger scalar("9");
    EC_Point res1, res2, res3;
    res1 = C.EC_Basic_Multiply( P1, scalar );
    res2 = C.EC_Double_and_Add( P1, scalar );
    res3 = C.EC_WindowNAF( P1, scalar, 4 );

    cout << endl;
    cout << "__________________________________________" << endl;
    cout << "Point P1 results : " << scalar <<" times " << P1 << endl;
    cout << "__________________________________________" << endl;
    cout << "Basic Multiply result : \t" << res1 << endl;
    cout << "Double and Add result : \t" << res2 << endl;
    cout << "WindowNAF result : \t\t" << res3 << endl;

    scalar = 5;
    res1 = C.EC_Basic_Multiply( P2, scalar );
    res2 = C.EC_Double_and_Add( P2, scalar );
    res3 = C.EC_WindowNAF( P2, scalar, 4 );

    cout << endl;
    cout << "__________________________________________" << endl;
    cout << "Point P2 results : " << scalar << " times " << P2 << endl;
    cout << "__________________________________________" << endl;
    cout << "Basic Multiply result : \t" << res1 << endl;
    cout << "Double and Add result : \t" << res2 << endl;
    cout << "WindowNAF result : \t\t" << res3 << endl;

    cout << endl;

    return 0;
}
