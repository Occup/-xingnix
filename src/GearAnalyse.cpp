#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <Eigen/Core>
#include <cmath>
#include <chrono>

using namespace std;
class frept
{
public:
  frept& operator=(const frept&);
  double fre() const {return _fre;}
  complex<double> Ri() const {return _Ri;}
  double A() const {return _A;}
  double P() const {return _P;}
  void construct(fstream&);
private:
  double _fre;
  double _A;
  double _P;
  complex<double> _Ri;
};
void frept::construct(fstream& dat)
{
  double A_in,W_in,t,y,r=0,i=0;
  dat >> A_in >> W_in;
  while(!dat.eof()){
    dat>>t>>y;
    r+=y*sin(W_in*t);
    i+=y*cos(W_in*t);
  }
  _fre=W_in;
  _Ri=(2*r/t,2*i/t);
  _A=abs(_Ri)/A_in;
  _P=arg(_Ri);
}
frept& frept::operator=(const frept& other)
{
  _fre=other.fre();
  _Ri=other.Ri();
  _A=other.A();
  _P=other.P();
  return *this;
}
// 曲线模型的顶点，模板参数：优化变量维度和数据类型
class CurveFittingVertex: public g2o::BaseVertex<10, Eigen::Matrix<double,1,10> > {
public:
  virtual void setToOriginImpl() // 重置
  {
    _estimate << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  }
  virtual void oplusImpl(const double* update) // 更新
  {
    _estimate += Eigen::Matrix<double,1,10>(update);
  }
  // 存盘和读盘：留空
  virtual bool read(istream& in) {}
  virtual bool write(ostream& out) const {}
};
// 误差模型 模板参数：观测值维度，类型，连接顶点类型
class CurveFittingEdge: public g2o::BaseUnaryEdge<1, double, CurveFittingVertex> {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW CurveFittingEdge(frept wp) :
    BaseUnaryEdge(), _pt(wp) {}
  // 计算曲线模型误差
  void computeError() {
    const CurveFittingVertex* v = static_cast<const CurveFittingVertex*>(_vertices[0]);
    const Eigen::Matrix<double,1,10> cof = v->estimate();
    complex<double> jw(0,_pt.fre());
    _error(0, 0) = abs(_pt.Ri()*(1.0+cof[5]*jw+cof[6]*pow(jw,2)+cof[7]*pow(jw,3)+cof[8]*pow(jw,4)+cof[8]*pow(jw,4)+cof[9]*pow(jw,5))
    		-(cof[0]+cof[1]*jw+cof[2]*pow(jw,2)+cof[3]*pow(jw,3)+cof[4]*pow(jw,4)));
  }
  
  virtual bool read(istream& in) {}
  virtual bool write(ostream& out) const {}
public:
  frept _pt;
};

int main(int argc, char** argv) {

    fstream ss;
    frept point;
    vector<frept> list;
    double tmp[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int i=1;i<argc;++i) {
      cout << "ate file "<<argv[i]<<endl;
      ss.open(argv[i], ios::in);
      point.construct(ss);
      ss.close();
      cout << point.A() << " " << point.fre() <<endl;
      list.push_back(point);
    }
    typedef g2o::BlockSolver< g2o::BlockSolverTraits<10,1>>  MyBlockSolver;
    typedef g2o::LinearSolverDense<MyBlockSolver::PoseMatrixType> MyLinearSolver;
    
  // setup the solver
    g2o::SparseOptimizer optimizer;
    optimizer.setVerbose(true);
    /*
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg
      (
       g2o::make_unique<MyBlockSolver>
       (
	g2o::make_unique<MyLinearSolver>()
	)
       );
    g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton
      (
       g2o::make_unique<MyBlockSolver>
       (
	g2o::make_unique<MyLinearSolver>()
	)
       );
    */
    g2o::OptimizationAlgorithmDogleg* solver = new g2o::OptimizationAlgorithmDogleg
       (
       g2o::make_unique<MyBlockSolver>
       (
	g2o::make_unique<MyLinearSolver>()
	)
       );
    optimizer.setAlgorithm(solver);
    // 往图中增加顶点
    CurveFittingVertex* v = new CurveFittingVertex();
    v->setEstimate(Eigen::Matrix<double,1,10>(tmp));
    v->setId(0);
    optimizer.addVertex(v);
    
    // 往图中增加边
    for(int i=0;i<list.size();++i){
      CurveFittingEdge* edge = new CurveFittingEdge(list[i]);
      edge->setId(i);
      edge->setVertex(0, v);                // 设置连接的顶点
      //	edge->setMeasurement(0);      // 观测数值
      //	edge->setInformation(Eigen::Matrix<double,1,1>::Identity()); // 信息矩阵：协方差矩阵之逆
      optimizer.addEdge(edge);
    }
    // 执行优化
    
    cout << "start optimization" << endl;
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    optimizer.initializeOptimization();
    optimizer.optimize(100);
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> time_used = chrono::duration_cast<
      chrono::duration<double>>(t2 - t1);
    cout << "solve time cost = " << time_used.count() << " seconds. " << endl;
    
	// 输出优化值
    Eigen::Matrix<double,1,10> cof = v->estimate();
    cout << "estimated model: " <<endl
	 << cof[0]<<"+"<<cof[1]<<"*s+"<<cof[2]<<"*s^2+"<<cof[3]<<"*s^3+"<<cof[4]<<"*s^4"<< endl
	 <<"------------------------------------------------------------------"<<endl
	 <<"1+"<< cof[5]<<"*s+"<<cof[6]<<"*s^2+"<<cof[7]<<"*s^3+"<<cof[8]<<"*s^4+"<<cof[9]<<"*s^5" << endl;

    return 0;
}
