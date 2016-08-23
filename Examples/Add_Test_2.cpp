#include<iostream>
#include "BigInteger.hpp"
#include "ECPO.hpp"
using namespace std;

int main() {

    EC_Curve C( 23, 1, 0 );
    C.EC_Validate();
    EC_Point points[24];
    points[0].zero = true;
    points[1].zero = false; points[1].x = 0; points[1].y = 0;
    points[2].zero = false; points[2].x = 1; points[2].y = 5;
    points[3].zero = false; points[3].x = 1; points[3].y = 18;
    points[4].zero = false; points[4].x = 9; points[4].y = 5;
    points[5].zero = false; points[5].x = 9; points[5].y = 18;
    points[6].zero = false; points[6].x = 11; points[6].y = 10;
    points[7].zero = false; points[7].x = 11; points[7].y = 13;
    points[8].zero = false; points[8].x = 13; points[8].y = 5;
    points[9].zero = false; points[9].x = 13; points[9].y = 18;
    points[10].zero = false; points[10].x = 15; points[10].y = 3;
    points[11].zero = false; points[11].x = 15; points[11].y = 20;
    points[12].zero = false; points[12].x = 16; points[12].y = 8;
    points[13].zero = false; points[13].x = 16; points[13].y = 15;
    points[14].zero = false; points[14].x = 17; points[14].y = 10;
    points[15].zero = false; points[15].x = 17; points[15].y = 13;
    points[16].zero = false; points[16].x = 18; points[16].y = 10;
    points[17].zero = false; points[17].x = 18; points[17].y = 13;
    points[18].zero = false; points[18].x = 19; points[18].y = 1;
    points[19].zero = false; points[19].x = 19; points[19].y = 22;
    points[20].zero = false; points[20].x = 20; points[20].y = 4;
    points[21].zero = false; points[21].x = 20; points[21].y = 19;
    points[22].zero = false; points[22].x = 21; points[22].y = 6;
    points[23].zero = false; points[23].x = 21; points[23].y = 17;

    for( int i = 0; i < 24; ++i ) {
        for( int j = 0; j < 24; ++j ) {
            EC_Point res;
            if( i != j ) res = C.EC_Add(points[i],points[j]);
            else res = C.EC_Double(points[i]);
            cout << res << "\t";
        }
        cout << endl;
    }

    return 0;
}
