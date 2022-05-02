#include <iomanip>
#include <string>
#include <random>
#include "raisim/RaisimServer.hpp"
#include "raisim/World.hpp"
#include "TApplication.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"

#include "kinematics.h"
#include "adaptive.h"

using namespace std;

void quaToRpy(raisim::Mat<3, 3> orientation, Vec3& RPY)
{
    Eigen::Quaterniond q(orientation.e());
    // roll (x-axis rotation)
    double sinr_cosp = +2.0 * (q.w() * q.x() + q.y() * q.z());
    double cosr_cosp = +1.0 - 2.0 * (q.x() * q.x() + q.y() * q.y());
    RPY[0] = atan2(sinr_cosp, cosr_cosp);

    // pitch (y-axis rotation)
    double sinp = +2.0 * (q.w() * q.y() - q.z() * q.x());
    if (fabs(sinp) >= 1)
    RPY[1] = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
    else
    RPY[1] = asin(sinp);

    // yaw (z-axis rotation)
    double siny_cosp = +2.0 * (q.w() * q.z() + q.x() * q.y());
    double cosy_cosp = +1.0 - 2.0 * (q.y() * q.y() + q.z() * q.z());
    RPY[2] = atan2(siny_cosp, cosy_cosp);
}

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  //创建世界、导入wk3
  raisim::World world;
  world.setTimeStep(0.001);   // 仿真环境中的更新间隔
  world.addGround();
  // world.setGravity(Eigen::VectorXd::Zero(3));
  auto wk3 = world.addArticulatedSystem(binaryPath.getDirectory() + "../rsc/wk3-mpcturn/Wukong3_rsm.urdf");
  wk3->setName("wk3");

  //打印基本信息
  cout.setf(ios::left);
  int posDim = wk3->getGeneralizedCoordinateDim();
  int dof = wk3->getDOF();
  double totalMass = wk3->getTotalMass();
  cout << "posDim: " << posDim << "\ndof: " << dof << "\ntotalMass: " << totalMass << endl;
  auto body = wk3->getBodyNames();
  auto mass = wk3->getMass();
  auto compositeMass = wk3->getCompositeMass();
  for (int i=0; i<body.size(); i++) cout  << setw(5) << i  
                                          << " name: " << setw(15) << body[i] 
                                          << " mass: " << setw(9) << mass[i] << endl;

  double upperMass;
  for (int i=0; i<body.size()-12; i++) upperMass += mass[i];
  cout << "upperMass: " << upperMass << endl;
  auto g = world.getGravity();
  cout << "g: " << g.e().transpose() << endl;

  //初始姿态
  Eigen::VectorXd posInit(posDim);
  Eigen::VectorXd pgain(dof), dgain(dof);
  Eigen::VectorXd posTarget(posDim), velTarget(dof), gf(dof);
  posInit <<  0,0, 0.90, 1,0,0,0,
              0., 
              0,0,0,0, 0,0,0,0,
              0,0,-0.549191,1.15799,0,-0.608327,
              0,0,-0.549191,1.15799,0,-0.608327;
  wk3->setGeneralizedCoordinate(posInit);
  wk3->setGeneralizedVelocity(Eigen::VectorXd::Zero(dof)); 
  wk3->setGeneralizedForce(Eigen::VectorXd::Zero(dof));

  /// launch raisim server
  //  运行前最好先打开unity可视化界面
  raisim::RaisimServer server(&world);
  server.launchServer();
  server.focusOn(wk3);
  std::this_thread::sleep_for(std::chrono::milliseconds(2000));

  //快速站稳
  raisim::Mat<3, 3> bodyRotation;
  Vec3 theta;
  double L1=0.35,L2=0.35;
  double thetala,thetara,thetalk,thetark;
  double A,B,C,D,Q,R,S,T,E,V,W;
  Eigen::MatrixXd Jac(4,3);
  Vec4 f;
  Vec3 fd,pd,dpd,p,dp;
  Mat3 Kp,Kd;
  pd << 0.8, 0, 0;
  dpd << 0, 0, 0;
  Kp.setZero();
  Kd.setZero();
  Kp.diagonal() << 10,10,10;
  Kd.diagonal() << 1,1,1;

  for (int i=0; i<2000; i++) {
    if(i%200==0){
      auto time = world.getWorldTime();
      auto pos = wk3->getGeneralizedCoordinate().e();
      auto vel = wk3->getGeneralizedVelocity().e();
      wk3->getBodyOrientation(0, bodyRotation);
      quaToRpy(bodyRotation, theta);

      p << pos(2), pos(0), theta(1);
      dp << vel(2), vel(0), vel(4);

      fd = Kp*(pd-p) + Kd*(dpd-dp);

      thetara = pos(21);
      thetala = pos(27);
      thetark = pos(19);
      thetalk = pos(25);

      A = -L1*cos(thetala) - L2*cos(thetala+thetalk);
      B = -L1*sin(thetala) - L2*sin(thetala+thetalk);
      C = -L1*cos(thetara) - L2*cos(thetara+thetark);
      D = -L1*sin(thetara) - L2*sin(thetara+thetark);
      Q = -L2*cos(thetala+thetalk);
      R = -L2*sin(thetala+thetalk);
      S = -L2*cos(thetara+thetark);
      T = -L2*sin(thetara+thetark);  
      E = C*B - A*D;
      V = Q*B - R*A;
      W = S*D - T*C;

      Jac <<  C*V/E,  D*V/E,  (-V-Q*D+R*C)/(2*E)-1/2,
              0,      0,      -0.5,
              -A*W/E, -B*W/E, (W+S*B-T*A)/(2*E)-1/2,
              0,      0,      -0.5;
      
      f = Jac*fd;

      gf(24) = f(0);
      gf(23) = f(1);
      gf(18) = f(2);
      gf(17) = f(3);

      // wk3->setGeneralizedForce(gf);

      cout << "worldtime: " << time << endl;
      cout << "pos: " << pos.transpose() << endl;
      cout << "vel: " << vel.transpose() << endl;
      cout << "p: " << p << endl;
      cout << "dp: " << dp << endl;
      cout << "Jac: " << Jac << endl;
      cout << "fd: " << fd << endl;
      cout << "f: " << f << endl;
    }

    this_thread::sleep_for(chrono::microseconds(900));
    server.integrateWorldThreadSafe();
  }

  server.killServer();
}
