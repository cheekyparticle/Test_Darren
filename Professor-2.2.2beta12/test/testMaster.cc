// Test for ParamPoints class
#include "Professor/ParamPoints.h"
#include "Professor/ProfMaster.h"
#include <iostream>

using namespace std;

int main(int argc, char** argv) {

  const vector<double> anchor1{0,0}, anchor2{0,1}, anchor3{0,2};
  const vector< vector<double> > points{anchor1, anchor2, anchor3};
  const vector<double> val{0,1,2};
  const Professor::ParamPoints pp(points);

  // Professor::ProfMaster pm(pp);
  // pm.addIpol("Test1", val, 1);
  // // pm.addIpol("Test2: 2 1 1.11022e-16 0 1");
  // for (size_t i = 0; i < 10; ++i) {
  //   double y = i*0.2;
  //   cout << "Value of Test1 at (x,y) = (0," << y << ") = " << pm.value("Test1", {0,y}) << endl;
  //   cout << "Value of Test2 at (x,y) = (0," << y << ") = " << pmvec.value("Test2", {0,y}) << endl;
  // }

  Professor::ProfMaster pmvec(points);
  pmvec.addIpol("Test1", val, 1);
  pmvec.addIpol("Test2: 2 1 1.11022e-16 0 1");
  for (size_t i = 0; i < 10; ++i) {
    double y = i*0.2;
    cout << "Value of Test1 at (x,y) = (0," << y << ") = " << pmvec.value("Test1", {0,y}) << endl;
    cout << "Value of Test2 at (x,y) = (0," << y << ") = " << pmvec.value("Test2", {0,y}) << endl;
  }

  return 0;
}
