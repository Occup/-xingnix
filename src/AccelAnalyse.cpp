#include <iostream>
#include <fstream>
//#include <sstream>
#include <ceres/ceres.h>

using namespace std;

struct CURVE_FITTING_COST
{
  CURVE_FITTING_COST(double Ax,double Ay,double Az,
		     double ax,double ay,double az):
    _Ax(Ax),_Ay(Ay),_Az(Az),
    _ax(ax),_ay(ay),_az(az)
  {}
  template <typename T>
  bool operator() (
		   const T * Mat, //3x5 Matrix
		   T * residual) const
  {
    residual[0] = T(_Ax) - (Mat[0]  + Mat[1]  * T(_ax) + Mat[2]  * T(_ay) + Mat[3]  * T(_az) + Mat[4]  * T(_ax) * T(_ax));
    residual[1] = T(_Ay) - (Mat[5]  + Mat[6]  * T(_ax) + Mat[7]  * T(_ay) + Mat[8]  * T(_az) + Mat[9]  * T(_ay) * T(_ay));
    residual[2] = T(_Az) - (Mat[10] + Mat[11] * T(_ax) + Mat[12] * T(_ay) + Mat[13] * T(_az) + Mat[14] * T(_az) * T(_az));
    return true;
  }
  const double _Ax,_Ay,_Az,_ax,_ay,_az;
};

int main(int args,char ** argv)
{
  /*
  stringstream ss;
  for (int i=1;i<args;++i) ss<<argv[i]<<" ";
  */
  double Ax,Ay,Az,ax,ay,az;
  double Mat[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  fstream data;
  if (args==2&&argv[1]!=NULL) data.open(argv[1], ios_base::in);
  else return 0;
  ceres::Problem problem;
  while (!data.eof()){
    data >>Ax>>Ay>>Az>>ax>>ay>>az;
    problem.AddResidualBlock(
			     new ceres::AutoDiffCostFunction<CURVE_FITTING_COST,3,15>(new CURVE_FITTING_COST(Ax,Ay,Az,ax,ay,az)),
			     NULL,
			     Mat);
  }
  data.close();
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;  // 增量方程如何求解
  options.minimizer_progress_to_stdout = true;   // 输出到cout
  ceres::Solver::Summary summary;                // 优化信息
  ceres::Solve(options,&problem,&summary);
  cout<<summary.BriefReport()<<endl;
  for (int i=0;i<3;++i) {
    for (int j=0;j<5;++j) cout<<Mat[5*i+j]<<"\t";
    cout<<endl;
  }
  return 0;
}
