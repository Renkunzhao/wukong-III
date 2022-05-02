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
  for (int i=0; i<body.size(); i++) cout  << setw(5) << i  
                                          << " name: " << setw(15) << body[i] 
                                          << " mass: " << setw(9) << mass[i] << endl;
  auto g = world.getGravity();
  Vec6 gTorso, fTorso, load;
  gTorso << (totalMass-mass[15]-mass[21])*g.e(),0,0,0;
  fTorso.setZero();
  load.setZero();
  cout << "gTorso: " << gTorso.transpose() << endl;

  //初始化WKLegKinematics类，根据需要的质心高度计算各个关节角度
  Eigen::VectorXd posTarget(posDim), velTarget(dof), feedForwardF(dof);
  posTarget <<  0,0, 0.75, 1,0,0,0,
                0., 
                0,0,0,0, 0,0,0,0,
                0., 0.042578, -0.492781, 0.97759, -0.03982, -0.535,
                0., -0.042567, -0.498771, 0.97759, 0.03964, -0.535;
  WKLegKinematics wkKin(posTarget);
  Vec3 endPr, endPl;
  Mat3 endR;
  // endPr << 0.023,-0.122,-0.8;
  // endPl << 0.023,0.122,-0.8;
  endPr << -0.1,-0.122,-0.7;
  endPl << 0.2,0.122,-0.7;
  endR.setIdentity();
  wkKin.inverseKin(endPr, endR, "RightLeg");
  wkKin.inverseKin(endPl, endR, "LeftLeg");
  for(int i=0;i<12;i++){
    posTarget[i+16] = wkKin.wkLink[i+1].q;
  }
  cout << "posTarget: " << posTarget.tail(12).transpose() << endl;
  velTarget.setZero();
  wk3->setGeneralizedCoordinate(posTarget);
  wk3->setGeneralizedVelocity(Eigen::VectorXd::Zero(dof)); 

  //控制参数
  Eigen::VectorXd pgain(dof), dgain(dof);
  wk3->setControlMode(raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE);
  pgain <<  0,0,0,0,0,0, 450,  80, 80, 80, 80, 80, 80, 80, 80, 
            1000, 1000, 200, 200, 1000, 200,
            1000, 1000, 200, 200, 1000, 200;
  dgain <<  0,0,0,0,0,0, 1.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
            1.5, 0.5, 25, 25, 0.3, 25,
            1.5, 0.5, 25, 25, 0.3, 25;
  // pgain <<  0,0,0,0,0,0, 450,  80, 80, 80, 80, 80, 80, 80, 80, 
  //           1000, 1000, 1000, 1000, 1000, 1000,
  //           1000, 1000, 1000, 1000, 1000, 1000;
  // dgain <<  0,0,0,0,0,0, 1.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
  //           25, 25, 25, 25, 25, 25,
  //           25, 25, 25, 25, 25, 25;
  wk3->setPdGains(pgain, dgain);
  wk3->setPdTarget(posTarget, velTarget);

  //初始化自适应控制参数
  WKAdaptive rAdapt, lAdapt;

  /// launch raisim server
  //  运行前最好先打开unity可视化界面
  raisim::RaisimServer server(&world);
  server.launchServer();
  server.focusOn(wk3);
  std::this_thread::sleep_for(std::chrono::milliseconds(2000));

  //调用ROOT环境准备画图
  TApplication app("app", &argc, argv);
  TCanvas* c = new TCanvas("c", "legTorque", 0, 0, 800, 600);
  TMultiGraph *mg1 = new TMultiGraph();
  auto gr1_1 = new TGraph();
  auto gr1_2 = new TGraph();
  auto gr1_3 = new TGraph();
  TMultiGraph *mg2 = new TMultiGraph();
  auto gr2_1 = new TGraph();
  auto gr2_2 = new TGraph();
  auto gr2_3 = new TGraph();
  TMultiGraph *mg3 = new TMultiGraph();
  auto gr3_1 = new TGraph();
  auto gr3_2 = new TGraph();
  auto gr3_3 = new TGraph();
  TMultiGraph *mg4 = new TMultiGraph();
  auto gr4_1 = new TGraph();
  auto gr4_2 = new TGraph();
  auto gr4_3 = new TGraph();
  TMultiGraph *mg5 = new TMultiGraph();
  auto gr5_1 = new TGraph();
  auto gr5_2 = new TGraph();
  auto gr5_3 = new TGraph();
  TMultiGraph *mg6 = new TMultiGraph();
  auto gr6_1 = new TGraph();
  auto gr6_2 = new TGraph();
  auto gr6_3 = new TGraph();

  TCanvas* c1 = new TCanvas("c1", "CoM", 0, 0, 1400, 600);
  TMultiGraph *mgCoM = new TMultiGraph();
  auto grCoMx = new TGraph();
  auto grCoMy = new TGraph();
  auto grCoMz = new TGraph();  
  auto grCoMtheta = new TGraph();  

  //快速站稳
  raisim::Mat<3, 3> bodyRotation;
  Vec3 thetad,theta;
  thetad << 0,0.2,0;
  for (int i=0; i<9000; i++) {
    auto pos = wk3->getGeneralizedCoordinate().e();
    auto vel = wk3->getGeneralizedVelocity().e();
    wk3->getBodyOrientation(0, bodyRotation);
    quaToRpy(bodyRotation, theta);
    wkKin.updateLinkq(pos);

    //VMC
    fTorso << 0,0,0,0,2500*(thetad[1]-theta[1]) + 20*(0-vel(4)),0;
    auto rtau = wkKin.RJac.transpose()*(-fTorso/2);
    auto ltau = wkKin.LJac.transpose()*(-fTorso/2);
    feedForwardF << Eigen::VectorXd::Zero(dof-12), rtau, ltau;
    // wk3->setGeneralizedForce(feedForwardF);

    //根据雅可比矩阵计算关节力矩
    fTorso << 0,0,100*(0.8-pos(2)) + 10*(0-vel(2)),
              1500*(thetad[0]-theta[0]) + 20*(0-vel(3)),
              2500*(thetad[1]-theta[1]) + 100*(0-vel(4)),
              0;
    wkKin.RTorque = wkKin.RJac.transpose()*((gTorso-fTorso+load)/2);
    wkKin.LTorque = wkKin.LJac.transpose()*((gTorso-fTorso+load)/2);
    feedForwardF << Eigen::VectorXd::Zero(dof-12), 1.05*wkKin.RTorque, 1.05*wkKin.LTorque;
    wk3->setGeneralizedForce(feedForwardF);

    //自适应控制算法
    rAdapt.updateq(pos.block(18,0,2,1), vel.block(17,0,2,1), posTarget.block(18,0,2,1), velTarget.block(17,0,2,1));
    lAdapt.updateq(pos.block(24,0,2,1), vel.block(23,0,2,1), posTarget.block(24,0,2,1), velTarget.block(23,0,2,1));
    fTorso << 0,0,0,0,500*(thetad[0]-theta[1]) + 20*(0-vel(4)),0;
    auto rAdaptTau = wkKin.RJac.transpose()*(-fTorso/2);
    auto lAdaptTau = wkKin.LJac.transpose()*(-fTorso/2);
    feedForwardF << Eigen::VectorXd::Zero(dof-12),  0,0,rAdapt.tau[0],rAdapt.tau[1],0,rAdaptTau[5], 
                                                    0,0,lAdapt.tau[0],lAdapt.tau[1],0,rAdaptTau[5];
    // wk3->setGeneralizedForce(feedForwardF);

    if(i%250==0){
      auto time = world.getWorldTime();
      auto pos = wk3->getGeneralizedCoordinate().e();
      auto vel = wk3->getGeneralizedVelocity().e();
      auto force = wk3->getGeneralizedForce().e();

      cout << "worldtime: " << time << endl;
      // cout << "posTarget: " << posTarget.tail(12).transpose() << endl;
      cout << "rTorque: " << force.block(15,0,6,1).transpose() << endl;
      cout << "lTorque: " << force.block(21,0,6,1).transpose() << endl;
      // cout << "rAdapt.tau: " << rAdapt.tau.transpose() << endl;
      // cout << "lAdapt.tau: " << lAdapt.tau.transpose() << endl;
      // cout << "Fd: " << Fd << endl;
      // cout << "Tdx: " << Tdx << endl;
      // cout << "Tdy: " << Tdy << endl;
      // cout << "Tdz: " << Tdz << endl;
      // cout << "wkKin.RTorque: " << wkKin.RTorque.transpose() << endl;
      // cout << "wkKin.LTorque: " << wkKin.LTorque.transpose() << endl;
      // cout << "r: right jacbian\n" << wkKin.RJac << endl;
      // cout << "r: right jacbian\n" << wkKin.LJac << endl;

      grCoMz->AddPoint(time, pos(2));
      grCoMtheta->AddPoint(time, theta(1));
      gr1_1->AddPoint(time, force(17));
      gr2_1->AddPoint(time, force(18));
      gr3_1->AddPoint(time, force(20));
      gr4_1->AddPoint(time, force(23));
      gr5_1->AddPoint(time, force(24));
      gr6_1->AddPoint(time, force(26));
      gr1_2->AddPoint(time, wkKin.RTorque(2));
      gr2_2->AddPoint(time, wkKin.RTorque(3));
      gr3_2->AddPoint(time, wkKin.RTorque(5));
      gr4_2->AddPoint(time, wkKin.LTorque(2));
      gr5_2->AddPoint(time, wkKin.LTorque(3));
      gr6_2->AddPoint(time, wkKin.LTorque(5));
      gr1_3->AddPoint(time, rAdapt.tau(0));
      gr2_3->AddPoint(time, rAdapt.tau(1));
      gr3_3->AddPoint(time, rAdaptTau[5]);
      gr4_3->AddPoint(time, lAdapt.tau(0));
      gr5_3->AddPoint(time, lAdapt.tau(1));
      gr6_3->AddPoint(time, lAdaptTau[5]);
    }

    // if(i>3000&&i<6000){
    //   load << 0,0,-300,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else{ 
    //   load << 0,0,0,0,0,0;
    //   wk3->clearExternalForcesAndTorques();
    // }

    // if(i>3000&&i<6000){
    //   load << 50,0,0,0,0,0;
    //   wk3->setExternalForce(1, load);
    // }
    // else{ 
    //   load << 0,0,0,0,0,0;
    //   wk3->clearExternalForcesAndTorques();
    // }

    if(i>3000&&i<6000){
      load << 0,0,-50,0,0,0;
      wk3->setExternalForce(13, load);
      wk3->setExternalForce(19, load);
    }
    else{ 
      load << 0,0,0,0,0,0;
      wk3->clearExternalForcesAndTorques();
    }

    // if(i>3000&&i<6000){
    //   load << 0,0,-50,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>6000&&i<9000){
    //   load << 0,0,-100,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>9000&&i<12000){
    //   load << 0,0,-150,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>12000&&i<15000){
    //   load << 0,0,-200,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>15000&&i<18000){
    //   load << 0,0,-250,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>18000&&i<21000){
    //   load << 0,0,-300,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else{
    //   load << 0,0,0,0,0,0;
    //   wk3->clearExternalForcesAndTorques();
    // }

    this_thread::sleep_for(chrono::microseconds(900));
    server.integrateWorldThreadSafe();
  }

  {//ROOT设置Pad等操作
    gr1_2->SetLineColor(2);
    gr2_2->SetLineColor(2);
    gr3_2->SetLineColor(2);
    gr4_2->SetLineColor(2);
    gr5_2->SetLineColor(2);
    gr6_2->SetLineColor(2);
    gr1_3->SetLineColor(4);
    gr2_3->SetLineColor(4);
    gr3_3->SetLineColor(4);
    gr4_3->SetLineColor(4);
    gr5_3->SetLineColor(4);
    gr6_3->SetLineColor(4);

    c1->Divide(2,1);
    c1->cd(1);
    grCoMz->Draw("AL");
    c1->cd(2);
    grCoMtheta->Draw("AL");

    c->Divide(3,2);
    c->cd(1);
    mg1->Add(gr1_1);
    mg1->Add(gr1_2);
    mg1->Add(gr1_3);
    mg1->Draw("AL");
    c->cd(2);
    mg2->Add(gr2_1);
    mg2->Add(gr2_2);
    mg2->Add(gr2_3);
    mg2->Draw("AL");
    c->cd(3);
    mg3->Add(gr3_1);
    mg3->Add(gr3_2);
    mg3->Add(gr3_3);
    mg3->Draw("AL");
    c->cd(4);
    mg4->Add(gr4_1);
    mg4->Add(gr4_2);
    mg4->Add(gr4_3);
    mg4->Draw("AL");
    c->cd(5);
    mg5->Add(gr5_1);
    mg5->Add(gr5_2);
    mg5->Add(gr5_3);
    mg5->Draw("AL");
    c->cd(6);
    mg6->Add(gr6_1);
    mg6->Add(gr6_2);
    mg6->Add(gr6_3);
    mg6->Draw("AL");

    app.Run();
  }

  server.killServer();
}
