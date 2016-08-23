#include<iostream>
#include "BigInteger.hpp"
#include "ECPO.hpp"
using namespace std;

int main() {

    EC_Curve C( 13, 3, 8 );
    C.EC_Validate();
    EC_Point points[9];
    points[0].zero = true;
    points[1].zero = false; points[1].x = 1; points[1].y = 5;
    points[2].zero = false; points[2].x = 1; points[2].y = 8;
    points[3].zero = false; points[3].x = 2; points[3].y = 3;
    points[4].zero = false; points[4].x = 2; points[4].y = 10;
    points[5].zero = false; points[5].x = 9; points[5].y = 6;
    points[6].zero = false; points[6].x = 9; points[6].y = 7;
    points[7].zero = false; points[7].x = 12; points[7].y = 2;
    points[8].zero = false; points[8].x = 12; points[8].y = 11;

    for( int i = 0; i < 9; ++i ) {
        for( int j = 0; j < 9; ++j ) {
            EC_Point res;
            if( i != j ) res = C.EC_Add(points[i],points[j]);
            else res = C.EC_Double(points[i]);
            cout << res << "\t";
        }
        cout << endl;
    }

    return 0;
}
