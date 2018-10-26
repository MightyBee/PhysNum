#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int main(){
	long double G(6.674E-11);
	long double mT(5.972E24);
	long double mL(7.342E22);
	long double z0(6378E3);
	long double zL(384400E3);
	long double zE(zL/(sqrt(mL/mT)+1));
	long double v0(sqrt((2*G)*(mT/z0 + mL/(zL-z0) - mT/zE - mL/(zL-zE))));
	cout << fixed << zE << endl << v0 << endl;
	return 0;
}
