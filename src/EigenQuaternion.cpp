#include <iostream>
#include <sstream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
using namespace std;
using namespace Eigen;

Quaterniond euler2Quaternion(const double Phi, const double Gamma, const double Psi)
{
    AngleAxisd PhiAngle(Phi, Vector3d::UnitZ());
    AngleAxisd PsiAngle(Psi, Vector3d::UnitY());
    AngleAxisd GammaAngle(Gamma, Vector3d::UnitX());

    Quaterniond q = PhiAngle * PsiAngle * GammaAngle;

    return q;
}
int main(int args,char ** argv)
{
  stringstream ss;
  double Phi,Psi,Gamma;
  Matrix3d R0,R1,R2;
  for (int i=1;i<args;++i) 
        ss<<argv[i]<<" ";
  cout <<ss.str()<<endl;
  cout <<ss.str().c_str()<<endl;
  ss>>Phi>>Psi>>Gamma;
  cout<<"EulerAngelRead:"
      <<Phi<<" "
      <<Psi<<" "
      <<Gamma<<endl;
  R0 << cos(Phi),  -sin(Phi),          0,
        sin(Phi),   cos(Phi),          0,
        0       ,   0       ,          1;
  R1 << cos(Psi),   0       ,   sin(Psi),
        0       ,   1       ,          0,
       -sin(Psi),   0       ,   cos(Psi);
  R2 << 1       ,   0       ,          0,
        0       , cos(Gamma),-sin(Gamma),
        0,        sin(Gamma), cos(Gamma);
  Quaterniond q = euler2Quaternion(Phi,Gamma,Psi);
  Matrix3d R = q.toRotationMatrix();
  Matrix3d r = R0*R1*R2;
  Vector3d V = R.eulerAngles(2,1,0);
  Vector3d v(atan2(r(1,0),r(0,0)),asin(-r(2,0)),atan2(r(2,1),r(2,2)));
  cout<<"QuaternionConverted:"
      <<q.coeffs().transpose()<<endl<<endl
      <<"from quaternion:"<<endl
      <<R<<endl<<endl
      <<"from eulerangle:"<<endl
      <<r<<endl<<endl
      <<"from quaternion:"<<endl
      <<V<<endl<<endl
      <<"from eulerangle:"<<endl
      <<v<<endl<<endl;
  return 0;
}
