#include<iostream>
#include<sys/time.h>
#include<cstdlib>
#include "BigInteger.hpp"
#include "ECPO.hpp"
using namespace std;

int main() {
    BigInteger p("416064700201658306196320137931"), A("3"), B("8");

// The curve and all the points have already been validated, so the statements have been removed or commented.

    EC_Curve C(p,A,B);
    // C.EC_Validate();

    EC_Point P1( BigInteger("5"), BigInteger("369763613166945291963257421414") );
    // cout << C.EC_Validate_Point( P1 ) << endl;
    EC_Point P2( BigInteger("9"), BigInteger("32949886211145258675155404324") );
    EC_Point P3( BigInteger("277589067115888363202510092850"), BigInteger("167916753052238204674165173582") );

    BigInteger scalar("123456789123456789");

    clock_t start1, end1, start2, end2;
    double cps = (double)(CLOCKS_PER_SEC);
    EC_Point res1, res2;

// Point P1

    // Double_and_Add measurement

        start1 = clock();
        res1 = C.EC_Double_and_Add( P1, scalar );
        end1 = clock();

    // WindowNAF measurement

        start2 = clock();
        res2 = C.EC_WindowNAF( P1, scalar, 4 );
        end2 = clock();

        cout << endl;
        cout << "__________________________________________" << endl;
        cout << "Point P1 results" << endl;
        cout << "__________________________________________" << endl;
        cout << "Double and Add running time : \t\t" << ( end1 - start1 ) / cps << endl;
        cout << "WindowNAF running time : \t\t" << ( end2 - start2 ) / cps << endl;
        cout << endl;
        cout << "Double and Add result : \t" << res1 << endl;
        cout << "WindowNAF result : \t\t" << res2 << endl;

// Point P2

    // Double_and_Add measurement

        start1 = clock();
        res1 = C.EC_Double_and_Add( P2, scalar );
        end1 = clock();

    // WindowNAF measurement

        start2 = clock();
        res2 = C.EC_WindowNAF( P2, scalar, 4 );
        end2 = clock();

        cout << endl << endl;
        cout << "__________________________________________" << endl;
        cout << "Point P2 results" << endl;
        cout << "__________________________________________" << endl;
        cout << "Double and Add running time : \t\t" << ( end1 - start1 ) / cps << endl;
        cout << "WindowNAF running time : \t\t" << ( end2 - start2 ) / cps << endl;
        cout << endl;
        cout << "Double and Add result : \t" << res1 << endl;
        cout << "WindowNAF result : \t\t" << res2 << endl;

// Point P3

    // Double_and_Add measurement

        start1 = clock();
        res1 = C.EC_Double_and_Add( P3, scalar );
        end1 = clock();

    // WindowNAF measurement

        start2 = clock();
        res2 = C.EC_WindowNAF( P3, scalar, 4 );
        end2 = clock();

        cout << endl << endl;
        cout << "__________________________________________" << endl;
        cout << "Point P3 results" << endl;
        cout << "__________________________________________" << endl;
        cout << "Double and Add running time : \t\t" << ( end1 - start1 ) / cps << endl;
        cout << "WindowNAF running time : \t\t" << ( end2 - start2 ) / cps << endl;
        cout << endl;
        cout << "Double and Add result : \t" << res1 << endl;
        cout << "WindowNAF result : \t\t" << res2 << endl;

    cout << endl;

    return 0;
}
