#include<iostream>
#include<utility>
#include<iomanip>
#include "BigInteger.hpp"
#include "ECPO.hpp"

#include<sys/time.h>
#include<cstdlib>

using namespace std;

int main() {
    BigInteger a("671998030559713968361666935769"),b("7867234521986094128367987634524"),c("786914378917592384751987982741");
    BigInteger d("327894769823456999348573946938457"), e("089762145337483191271827434177"), f("282174488599599500573849980909");
    BigInteger g("521419622856657689423872613771"), h("23411040409478237411787827481"), i("362736035870515331128527330659");
    cout << a << "\t" << BI_Miller_Rabin(a) << endl;
    cout << b << "\t" << BI_Miller_Rabin(b) << endl;
    cout << c << "\t" << BI_Miller_Rabin(c) << endl;
    cout << d << "\t" << BI_Miller_Rabin(d) << endl;
    cout << e << "\t" << BI_Miller_Rabin(e) << endl;
    cout << f << "\t" << BI_Miller_Rabin(f) << endl;
    cout << g << "\t" << BI_Miller_Rabin(g) << endl;
    cout << h << "\t" << BI_Miller_Rabin(h) << endl;
    cout << i << "\t" << BI_Miller_Rabin(i) << endl;
    return 0;
}
