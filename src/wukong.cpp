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
  Vec6 gTorso;
  gTorso << (totalMass-mass[15]-mass[21])*g.e(),0,0,0;
  cout << "gTorso: " << gTorso.transpose() << endl;

  auto massMatrix = wk3->getMassMatrix().e();
  cout << "MassMatrix:\n" << massMatrix << endl;

  //初始姿态
  Eigen::VectorXd posInit(posDim);
  Eigen::VectorXd pgain(dof), dgain(dof);
  Eigen::VectorXd posTarget(posDim), velTarget(dof), feedForwardF(dof);
  posInit <<  0,0, 0.90, 1,0,0,0,
              0., 
              0,0,0,0, 0,0,0,0,
              0., 0.042578, -0.492781, 0.97759, -0.03982, -0.535,
              0., -0.042567, -0.498771, 0.97759, 0.03964, -0.535;
  wk3->setGeneralizedCoordinate(posInit);
  wk3->setGeneralizedVelocity(Eigen::VectorXd::Zero(dof)); 
  // wk3->setIntegrationScheme(raisim::ArticulatedSystem::IntegrationScheme::SEMI_IMPLICIT);

  //初始化WKLegKinematics类，用于正逆运动学、雅可比矩阵计算
  WKLegKinematics wkKin(posInit);
  cout << "pTarget: " << wkKin.wkLink[12].p.transpose() << endl;                  
  //根据需要的质心高度计算各个关节角度
  Vec3 endPr, endPl;
  Mat3 endR;
  // endPr << 0.023,-0.122,-0.8;
  // endPl << 0.023,0.122,-0.8;
  // endPr << 0.0,-0.122,-0.8;
  // endPl << 0.0,0.122,-0.8;
  endPr << 0.0,-0.1,-0.8;
  endPl << 0.0,0.1,-0.8;
  // endPr << 0.2,-0.1,-0.5;
  // endPl << -0.2,0.1,-0.5;  
  endR.setIdentity();
  wkKin.inverseKin(endPr, endR, "RightLeg");
  wkKin.inverseKin(endPl, endR, "LeftLeg");
  //利用雅可比矩阵,根据机器人重量及负载计算关节力矩
  Vec6 load;
  load.setZero();

  //初始化自适应控制参数
  WKAdaptive rAdapt, lAdapt;

  //控制参数
  wk3->setControlMode(raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE);
  // wk3->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);
  // pgain <<  0,0,0,0,0,0, 450,  80, 80, 80, 80, 80, 80, 80, 80, 
  //           350, 500, 1200, 500, 100, 300,
  //           350, 500, 1200, 500, 100, 300;
  // dgain <<  0,0,0,0,0,0, 1.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
  //           1.5, 0.5, 45, 25, 0.3, 3,
  //           1.5, 0.5, 45, 25, 0.3, 3;
  // pgain <<  0,0,0,0,0,0, 450,  80, 80, 80, 80, 80, 80, 80, 80, 
  //           1000, 1000, 1000, 100, 1000, 1000,
  //           1000, 1000, 1000, 100, 1000, 1000;
  // dgain <<  0,0,0,0,0,0, 1.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
  //           1.5, 0.5, 45, 25, 0.3, 45,
  //           1.5, 0.5, 45, 25, 0.3, 45;
  pgain <<  0,0,0,0,0,0, 450,  80, 80, 80, 80, 80, 80, 80, 80, 
            1000, 1000, 200, 200, 1000, 200,
            1000, 1000, 200, 200, 1000, 200;
  dgain <<  0,0,0,0,0,0, 1.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
            1.5, 0.5, 25, 25, 0.3, 25,
            1.5, 0.5, 25, 25, 0.3, 25;
  posTarget = posInit;
  for(int i=0;i<12;i++){
    posTarget[i+16] = wkKin.wkLink[i+1].q;
  }
  cout << "posTarget: " << posTarget.tail(12).transpose() << endl;
  velTarget.setZero();
  feedForwardF << Eigen::VectorXd::Zero(dof-12), wkKin.RTorque, wkKin.LTorque;
  wk3->setPdGains(pgain, dgain);
  wk3->setPdTarget(posTarget, velTarget);
  // wk3->setGeneralizedForce(feedForwardF);

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

  TCanvas* c1 = new TCanvas("c1", "CoM", 0, 0, 800, 600);
  // c1->Divide(1,3);
  TMultiGraph *mgCoM = new TMultiGraph();
  auto grCoMx = new TGraph();
  auto grCoMy = new TGraph();
  auto grCoMz = new TGraph();  

  //快速站稳
  raisim::Mat<3, 3> bodyRotation;
  Vec3 eulAng;
  double HipPGain = -10.0;
  for (int i=0; i<21000; i++) {

    wk3->getBodyOrientation(0, bodyRotation);
    quaToRpy(bodyRotation, eulAng);
    posTarget[18] = posInit[18] + HipPGain*(0-eulAng[1]);
    posTarget[24] = posInit[24] + HipPGain*(0-eulAng[1]);//加个增益使pitch更稳
    posTarget[21] = posInit[21] + HipPGain*(0-eulAng[1]);
    posTarget[27] = posInit[27] + HipPGain*(0-eulAng[1]);//加个增益使pitch更稳
    wk3->setPdTarget(posTarget, velTarget);

    auto pos = wk3->getGeneralizedCoordinate().e();
    auto vel = wk3->getGeneralizedVelocity().e();

    //根据雅可比矩阵计算关节力矩
    wkKin.updateLinkq(pos);
    wkKin.RTorque = wkKin.RJac.transpose()*((gTorso+load)/2);
    wkKin.LTorque = wkKin.LJac.transpose()*((gTorso+load)/2);
    feedForwardF << Eigen::VectorXd::Zero(dof-12), 1.05*wkKin.RTorque, 1.05*wkKin.LTorque;
    // wk3->setGeneralizedForce(feedForwardF);

    //自适应控制算法
    rAdapt.updateq(pos.block(18,0,2,1), vel.block(17,0,2,1), posTarget.block(18,0,2,1), velTarget.block(17,0,2,1));
    lAdapt.updateq(pos.block(24,0,2,1), vel.block(23,0,2,1), posTarget.block(24,0,2,1), velTarget.block(23,0,2,1));
    feedForwardF << Eigen::VectorXd::Zero(dof-12),  0,0,rAdapt.tau[0],rAdapt.tau[1],0,0, 
                                                    0,0,lAdapt.tau[0],lAdapt.tau[1],0,0;
    wk3->setGeneralizedForce(feedForwardF);

    if(i%250==0){
      auto time = world.getWorldTime();
      auto pos = wk3->getGeneralizedCoordinate().e();
      auto vel = wk3->getGeneralizedVelocity().e();
      auto force = wk3->getGeneralizedForce().e();

      // wkKin.updateLinkq(pos);
      // //根据雅可比矩阵计算关节力矩
      // wkKin.jointTorque(gHipy, gTorso, load);
      // feedForwardF << Eigen::VectorXd::Zero(dof-12), wkKin.RTorque, wkKin.LTorque;
      // wk3->setGeneralizedForce(feedForwardF);

      cout << "worldtime: " << time << endl;
      cout << "posTarget: " << posTarget.tail(12).transpose() << endl;
      cout << "rAdapt.tau: " << force.block(15,0,6,1).transpose() << endl;
      cout << "lAdapt.tau: " << force.block(21,0,6,1).transpose() << endl;
      cout << "rAdapt.tau: " << rAdapt.tau.transpose() << endl;
      cout << "lAdapt.tau: " << lAdapt.tau.transpose() << endl;
      cout << "wkKin.RTorque: " << wkKin.RTorque.transpose() << endl;
      cout << "wkKin.LTorque: " << wkKin.LTorque.transpose() << endl;

      grCoMz->AddPoint(time, pos(2));
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
      // gr3_3->AddPoint(time, rAdapt.tau(2));
      gr4_3->AddPoint(time, lAdapt.tau(0));
      gr5_3->AddPoint(time, lAdapt.tau(1));
      // gr6_3->AddPoint(time, lAdapt.tau(2));
    }

    // if(i>3000&&i<6000){
    //   load << 0,0,-250,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else{ 
    //   load << 0,0,0,0,0,0;
    //   // wk3->setExternalForce(0, load);
    //   wk3->clearExternalForcesAndTorques();
    // }

    if(i>3000&&i<6000){
      load << 0,0,-50,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else if(i>6000&&i<9000){
      load << 0,0,-100,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else if(i>9000&&i<12000){
      load << 0,0,-150,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else if(i>12000&&i<15000){
      load << 0,0,-200,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else if(i>15000&&i<18000){
      load << 0,0,-250,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else if(i>18000&&i<21000){
      load << 0,0,-300,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else{
      load << 0,0,0,0,0,0;
      wk3->clearExternalForcesAndTorques();
    }

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

    grCoMz->Draw("AL");

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
