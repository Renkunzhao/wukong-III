#include <iomanip>
#include <string>
#include <random>
#include <vector>
#include <algorithm> 
#include "raisim/RaisimServer.hpp"
#include "raisim/World.hpp"
#include "matplotlibcpp.h"

#include "kinematics.h"
#include "adaptive3.h"

using namespace std;
namespace plt = matplotlibcpp;

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
  for (int i=0; i<body.size(); i++) cout  << setw(5) << i  
                                          << " name: " << setw(15) << body[i] 
                                          << " mass: " << setw(9) << mass[i] << endl;
  auto g = world.getGravity();
  Vec6 gTorso, fTorso, load, footFix;
  gTorso << (totalMass-mass[15]-mass[21])*g.e(),0,0,0;
  fTorso.setZero();
  load.setZero();
  cout << "gTorso: " << gTorso.transpose() << endl;

  //初始化WKLegKinematics类，根据需要的质心高度计算各个关节角度
  Eigen::VectorXd posInit(posDim), posTarget(posDim), velTarget(dof), feedForwardF(dof);
  Eigen::VectorXd pos(posDim), vel(dof), velLast(dof), force(dof);
  double CoMTheta = 0.1;
  double CoMz = 0.6;
  posTarget <<  0,0, CoMz+0.05, 1,0,0,0,
                0., 
                0,0,0,0, 0,0,0,0,
                0., 0.042578, -0.492781, 0.97759, -0.03982, -0.535,
                0., -0.042567, -0.498771, 0.97759, 0.03964, -0.535;
  WKLegKinematics wkKin(posTarget);
  Eigen::AngleAxisd pitchAngleCoM(CoMTheta,Eigen::Vector3d::UnitY());
  Eigen::Quaterniond quaternion(pitchAngleCoM);
  Eigen::AngleAxisd pitchAngle(-CoMTheta,Eigen::Vector3d::UnitY());
  Vec3 endPr, endPl;
  Mat3 endR;
  // endPr << 0.023,-0.122,-0.8;
  // endPl << 0.023,0.122,-0.8;
  endPr << -0.,-0.122,-CoMz;
  endPl << 0.,0.122,-CoMz;
  endPr = pitchAngle.matrix()*endPr;
  endPl = pitchAngle.matrix()*endPl;
  endR = pitchAngle.matrix();
  wkKin.inverseKin(endPr, endR, "RightLeg");
  wkKin.inverseKin(endPl, endR, "LeftLeg");
  for(int i=0;i<12;i++){
    posTarget[i+16] = wkKin.wkLink[i+1].q;
  }
  posTarget[3] = quaternion.coeffs()[3];
  posTarget[4] = quaternion.coeffs()[0];
  posTarget[5] = quaternion.coeffs()[1];
  posTarget[6] = quaternion.coeffs()[2];
  posInit = posTarget;
  cout << "posTarget: " << posTarget.tail(12).transpose() << endl;
  velTarget.setZero();
  wk3->setGeneralizedCoordinate(posTarget);
  wk3->setGeneralizedVelocity(Eigen::VectorXd::Zero(dof)); 

  //控制参数
  Eigen::VectorXd pgain(dof), dgain(dof);
  wk3->setControlMode(raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE);
  pgain <<  0,0,0,0,0,0, 450,  80, 80, 80, 80, 80, 80, 80, 80, 
            1000, 1000, 400, 200, 1000, 400,
            1000, 1000, 400, 200, 1000, 400;
  dgain <<  0,0,0,0,0,0, 1.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
            1.5, 0.5, 40, 25, 0.3, 80,
            1.5, 0.5, 40, 25, 0.3, 80;
  // pgain <<  0,0,0,0,0,0, 450,  80, 80, 80, 80, 80, 80, 80, 80, 
  //           350, 500, 1200, 500, 100, 300,
  //           350, 500, 1200, 500, 100, 300;
  // dgain <<  0,0,0,0,0,0, 1.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
  //           1.5, 0.5, 45, 25, 0.3, 3,
  //           1.5, 0.5, 45, 25, 0.3, 3;
  wk3->setPdGains(pgain, dgain);
  wk3->setPdTarget(posTarget, velTarget);

  //初始化自适应控制参数
  WKAdaptive3 rAdapt, lAdapt;

  /// launch raisim server
  //  运行前最好先打开unity可视化界面
  raisim::RaisimServer server(&world);
  server.launchServer();
  server.focusOn(wk3);
  std::this_thread::sleep_for(std::chrono::milliseconds(2000));

  //快速站稳
  raisim::Mat<3, 3> bodyRotation;
  Vec3 thetad,theta;
  thetad << 0,CoMTheta,0;
  for (int i=0; i<2000; i++) {
    auto pos = wk3->getGeneralizedCoordinate().e();
    auto vel = wk3->getGeneralizedVelocity().e();
    wk3->getBodyOrientation(0, bodyRotation);
    quaToRpy(bodyRotation, theta);
    wkKin.updateLinkq(pos);

    //根据雅可比矩阵计算关节力矩
    fTorso << 0,0,0*(CoMz+0.045-pos(2)) + 0*(0-vel(2)),
              0*(thetad[0]-theta[0]) + 0*(0-vel(3)),
              2500*(thetad[1]-theta[1]) + 20*(0-vel(4)),
              0;
    wkKin.RTorque = wkKin.RJac.transpose()*((gTorso-fTorso+load)/2);
    wkKin.LTorque = wkKin.LJac.transpose()*((gTorso-fTorso+load)/2);
    feedForwardF << Eigen::VectorXd::Zero(dof-12), wkKin.RTorque, wkKin.LTorque;
    wk3->setGeneralizedForce(feedForwardF);

    this_thread::sleep_for(chrono::microseconds(900));
    server.integrateWorldThreadSafe();
  }

  pos.setZero();
  vel.setZero();
  velLast.setZero();
  force.setZero();
  vector<double> pltt,pltz,plttheta,q1,q2,q3,a1,a2,a3,a4,a5,a6,a7,a8,a9,f1,f2,f3,fJac1,fJac2,fJac3,fAda1,fAda2,fAda3;
  //增加负载
  for (int i=0; i<21000; i++) {
    auto time = world.getWorldTime();
    pos = wk3->getGeneralizedCoordinate().e();
    velLast = vel;
    vel = wk3->getGeneralizedVelocity().e();
    force = wk3->getGeneralizedForce().e();
    wk3->getBodyOrientation(0, bodyRotation);
    quaToRpy(bodyRotation, theta);
    wkKin.updateLinkq(pos);

    //锁住脚底
    // footFix << 0,0,-1000,0,0,0;
    // wk3->setExternalForce(15,footFix);
    // wk3->setExternalForce(21,footFix);

    //抑制晃动
    // posTarget[18] = posInit[18] - 10*(0-theta[1]);
    // posTarget[24] = posInit[24] - 10*(0-theta[1]);//加个增益使pitch更稳
    // posTarget[21] = posInit[21] - 10*(0-theta[1]);
    // posTarget[27] = posInit[27] - 10*(0-theta[1]);//加个增益使pitch更稳
    // wk3->setPdTarget(posTarget, velTarget);

    //VMC+Jac 无负载前馈
    // fTorso << 0,0,0,0,2500*(thetad[1]-theta[1]) + 20*(0-vel(4)),0;
    // auto rtau = wkKin.RJac.transpose()*(-fTorso/2);
    // auto ltau = wkKin.LJac.transpose()*(-fTorso/2);
    // feedForwardF << Eigen::VectorXd::Zero(dof-12), rtau, ltau;
    // wk3->setGeneralizedForce(feedForwardF);

    //VMC+Jac 有负载前馈
    // fTorso << 0,0,0*(CoMz+0.045-pos(2)) + 0*(0-vel(2)),
    //           0*(thetad[0]-theta[0]) + 0*(0-vel(3)),
    //           2500*(thetad[1]-theta[1]) + 20*(0-vel(4)),
    //           0;
    // wkKin.RTorque = wkKin.RJac.transpose()*((gTorso-fTorso+load)/2);
    // wkKin.LTorque = wkKin.LJac.transpose()*((gTorso-fTorso+load)/2);
    // feedForwardF << Eigen::VectorXd::Zero(dof-12), wkKin.RTorque, wkKin.LTorque;
    // wk3->setGeneralizedForce(feedForwardF);

    //Adaptive
    if((vel.tail(12)-velLast.tail(12)).array().abs().maxCoeff()<0.1){
      rAdapt.q << pos(21), pos(19), pos(18);
      rAdapt.qd << posTarget(21), posTarget(19), posTarget(18);
      rAdapt.dq << vel(20), vel(18), vel(17);
      rAdapt.dqd << velTarget(20), velTarget(18), velTarget(17);
      lAdapt.q << pos(27), pos(25), pos(24);
      lAdapt.qd << posTarget(27), posTarget(25), posTarget(24);
      lAdapt.dq << vel(26), vel(24), vel(23);
      lAdapt.dqd << velTarget(26), velTarget(24), velTarget(23);
      rAdapt.updates();
      lAdapt.updates();
      feedForwardF << Eigen::VectorXd::Zero(dof-12),  0,0,rAdapt.tau[2],rAdapt.tau[1],0,rAdapt.tau[0], 
                                                      0,0,lAdapt.tau[2],lAdapt.tau[1],0,lAdapt.tau[0];
      wk3->setGeneralizedForce(feedForwardF);
    }

    if(i%100==0){
      cout << "worldtime: " << time << endl;
      // cout << "rAdapt.Y:\n" << rAdapt.Y << endl;
      // cout << "lAdapt.Y:\n" << lAdapt.Y << endl;

      pltt.push_back(time);
      pltz.push_back(pos(2));
      plttheta.push_back(theta(1));
      q1.push_back(pos(27));
      q2.push_back(pos(25));
      q3.push_back(pos(24));
      a1.push_back(lAdapt.a_(0));
      a2.push_back(lAdapt.a_(1));
      a3.push_back(lAdapt.a_(2));
      a4.push_back(lAdapt.a_(3));
      a5.push_back(lAdapt.a_(4));
      a6.push_back(lAdapt.a_(5));
      a7.push_back(lAdapt.a_(6));
      a8.push_back(lAdapt.a_(7));
      a9.push_back(lAdapt.a_(8));
      // q1.push_back(pos(21));
      // q2.push_back(pos(19));
      // q3.push_back(pos(18));
      // a1.push_back(lAdapt.a_(0));
      // a2.push_back(lAdapt.a_(1));
      // a3.push_back(lAdapt.a_(2));
      // a4.push_back(lAdapt.a_(3));
      // a5.push_back(lAdapt.a_(4));
      // a6.push_back(lAdapt.a_(5));
      // a7.push_back(lAdapt.a_(6));
      // a8.push_back(lAdapt.a_(7));
      // a9.push_back(lAdapt.a_(8));
      f1.push_back(force(20));
      f2.push_back(force(18));
      f3.push_back(force(17));
      fJac1.push_back(wkKin.RTorque(5));
      fJac2.push_back(wkKin.RTorque(3));
      fJac3.push_back(wkKin.RTorque(2));
      fAda1.push_back(rAdapt.tau(0));
      fAda2.push_back(rAdapt.tau(1));
      fAda3.push_back(rAdapt.tau(2));
    }

    if(i==1900||i==9900){
      cout << "CoMz: " << pltz.back() << endl;
      cout << "CoMtheta: " << plttheta.back() << endl;
      cout << "HipY: " << q3.back() << endl;
      cout << "KneeY: " << q2.back() << endl;
      cout << "AnkleY: " << q1.back() << endl;
    }

    // if(i>2000&&i<10000){
    //   load << 0,0,-300,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else{ 
    //   load << 0,0,0,0,0,0;
    //   wk3->clearExternalForcesAndTorques();
    // }

    // if(i>3000&&i<6000){
    //   load << -50,0,0,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>6000&&i<9000){
    //   load << -100,0,0,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else{
    //   load << 0,0,0,0,0,0;
    //   wk3->clearExternalForcesAndTorques();
    // }

    if(i>3000&&i<6000){
      load << 0,0,-50,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else if(i>6000&&i<9000){
      load << 0,0,-150,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else if(i>9000&&i<12000){
      load << 0,0,-250,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else if(i>12000&&i<15000){
      load << 0,0,-350,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else if(i>15000&&i<18000){
      load << 0,0,-450,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else if(i>18000&&i<21000){
      load << 0,0,-550,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else{
      load << 0,0,0,0,0,0;
      wk3->clearExternalForcesAndTorques();
    }

    this_thread::sleep_for(chrono::microseconds(900));
    server.integrateWorldThreadSafe();
  }

  vector<double> pltzd(pltt.size(), CoMz+0.045);
  vector<double> pltthetad(pltt.size(), CoMTheta);
  vector<double> q1d(pltt.size(), posTarget(21));
  vector<double> q2d(pltt.size(), posTarget(19));
  vector<double> q3d(pltt.size(), posTarget(18));

  cout << "CoMzd: " << pltzd.back() << endl;
  cout << "CoMthetad: " << pltthetad.back() << endl;
  cout << "HipYd: " << q3d.back() << endl;
  cout << "KneeYd: " << q2d.back() << endl;
  cout << "AnkleYd: " << q1d.back() << endl;

  plt::figure(1);
  plt::plot(pltt, pltz, "r-",{{"label", "z"}});
  plt::plot(pltt, pltzd, "k:",{{"label", "z$_{set}$"}});
  plt::legend();
  // plt::ylim(0.7, 0.8);
  plt::xlabel("时间(s)");
  plt::ylabel("高度(m)");
  plt::title("机器人质心高度$z$");

  plt::figure(2);
  plt::plot(pltt, plttheta, "r-", {{"label", "$\\theta$"}});
  plt::plot(pltt, pltthetad, "k:", {{"label", "$\\theta_{set}$"}});
  plt::legend();
  plt::xlabel("时间(s)");
  plt::ylabel("角度(rad)");
  plt::title("机器人俯仰角度$\\theta$");

  plt::figure(3);
  plt::plot(pltt, q1, "r-", {{"label", "AnkleY"}});
  plt::plot(pltt, q1d, "k:", {{"label", "AnkleY$_{set}$"}});
  plt::legend();
  plt::xlabel("时间(s)");
  plt::ylabel("角度(rad)");
  plt::title("关节角度AnkleY");

  plt::figure(4);
  plt::plot(pltt, q2, "g-", {{"label", "KneeY"}});
  plt::plot(pltt, q2d, "k:", {{"label", "KneeY$_{set}$"}});
  plt::legend();
  plt::xlabel("时间(s)");
  plt::ylabel("角度(rad)");
  plt::title("关节角度KneeY");

  plt::figure(5);
  plt::plot(pltt, q3, "m-", {{"label", "HipY"}});
  plt::plot(pltt, q3d, "k:", {{"label", "HipY$_{set}$"}});
  plt::legend();
  plt::xlabel("时间(s)");
  plt::ylabel("角度(rad)");
  plt::title("关节角度HipY");

  plt::figure(6);
  plt::plot(pltt, a1, "b-", {{"label", "a1"}});
  plt::plot(pltt, a2, "c-", {{"label", "a2"}});
  plt::plot(pltt, a3, "k-", {{"label", "a3"}});
  plt::plot(pltt, a4, "g-", {{"label", "a4"}});
  plt::plot(pltt, a5, "m-", {{"label", "a5"}});
  plt::plot(pltt, a6, "r-", {{"label", "a6"}});
  plt::legend();
  plt::xlabel("时间(s)");
  plt::title("动力学参数$a$");

  plt::figure(7);
  plt::plot(pltt, a7, "k-", {{"label", "a7"}});
  plt::plot(pltt, a8, "g-", {{"label", "a8"}});
  plt::plot(pltt, a9, "m-", {{"label", "a9"}});
  plt::legend();
  plt::xlabel("时间(s)");
  plt::title("动力学参数$a$");

  plt::figure(8);
  plt::plot(pltt, f3, "k:", {{"label", "f$_{real}$"}});
  // plt::plot(pltt, fJac3, "r-", {{"label", "f$_{Jac}$"}});
  plt::plot(pltt, fAda3, "r-", {{"label", "f$_{Ada}$"}});
  plt::legend();
  plt::xlabel("时间(s)");
  plt::ylabel("力矩($N\\cdot m$)");
  plt::title("HipY关节力矩$\\tau$");

  plt::figure(9);
  plt::plot(pltt, f2, "k:", {{"label", "f$_{real}$"}});
  // plt::plot(pltt, fJac2, "r-", {{"label", "f$_{Jac}$"}});
  plt::plot(pltt, fAda2, "r-", {{"label", "f$_{Ada}$"}});
  plt::legend();
  plt::xlabel("时间(s)");
  plt::ylabel("力矩($N\\cdot m$)");
  plt::title("KneeY关节力矩$\\tau$");

  plt::figure(10);
  plt::plot(pltt, f1, "k:", {{"label", "f$_{real}$"}});
  // plt::plot(pltt, fJac1, "r-", {{"label", "f$_{Jac}$"}});
  plt::plot(pltt, fAda1, "r-", {{"label", "f$_{Ada}$"}});
  plt::legend();
  plt::xlabel("时间(s)");
  plt::ylabel("力矩($N\\cdot m$)");
  plt::title("AnkleY关节力矩$\\tau$");

  plt::show();

  server.killServer();
}
