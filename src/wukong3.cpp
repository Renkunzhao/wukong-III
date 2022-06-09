#include <iomanip>
#include <string>
#include <random>
#include <vector>
#include <algorithm> 
#include "raisim/RaisimServer.hpp"
#include "raisim/World.hpp"
#include "TApplication.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"

#include "kinematics.h"
#include "adaptive3.h"

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

double filter(vector<double>& q){
  for(int i=0;i<q.size();i++){
    cout << q[i] << " ";
  }
  cout << endl;
  if(q.size()>15){
    sort(q.begin(), q.end());
    double sum = accumulate(q.begin()+5, q.end()-5, 0.0);
    double mean =  sum / (q.size()-10);
    q.clear();
    cout << mean << endl;
    return mean;
  }
  else{
    cout << q.size() << endl;
    return 0;
  }
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
  double CoMz = 0.7;
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
            1000, 1000, 200, 200, 1000, 200,
            1000, 1000, 200, 200, 1000, 200;
  dgain <<  0,0,0,0,0,0, 1.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
            1.5, 0.5, 25, 25, 0.3, 25,
            1.5, 0.5, 25, 25, 0.3, 25;
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

  TCanvas* c2 = new TCanvas("c2", "dynamics", 0, 0, 1400, 600);
  TMultiGraph *mgA1 = new TMultiGraph();
  auto gra1 = new TGraph();
  auto gra2 = new TGraph();
  auto gra3 = new TGraph();
  auto gra4 = new TGraph();
  auto gra5 = new TGraph();
  auto gra6 = new TGraph();
  TMultiGraph *mgA2 = new TMultiGraph();
  auto gra7 = new TGraph();
  auto gra8 = new TGraph();
  auto gra9 = new TGraph();

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
    fTorso << 0,0,0*(CoMz+0.044-pos(2)) + 0*(0-vel(2)),
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
  //增加负载
  for (int i=0; i<1000; i++) {
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
    // fTorso << 0,0,0*(CoMz+0.044-pos(2)) + 0*(0-vel(2)),
    //           0*(thetad[0]-theta[0]) + 0*(0-vel(3)),
    //           2500*(thetad[1]-theta[1]) + 20*(0-vel(4)),
    //           0;
    // wkKin.RTorque = wkKin.RJac.transpose()*((-fTorso+load)/2);
    // wkKin.LTorque = wkKin.LJac.transpose()*((-fTorso+load)/2);
    // feedForwardF << Eigen::VectorXd::Zero(dof-12), wkKin.RTorque, wkKin.LTorque;
    // wk3->setGeneralizedForce(feedForwardF);

    if((vel.tail(12)-velLast.tail(12)).array().abs().maxCoeff()<0.1){
      //自适应控制算法
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
      cout << "rAdapt.a: " << rAdapt.a_.transpose() << endl;
      cout << "lAdapt.a: " << lAdapt.a_.transpose() << endl;
      cout << rAdapt.tau.transpose() << endl;

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
      gr1_3->AddPoint(time, rAdapt.tau(2));
      gr2_3->AddPoint(time, rAdapt.tau(1));
      gr3_3->AddPoint(time, rAdapt.tau(0));
      gr4_3->AddPoint(time, lAdapt.tau(2));
      gr5_3->AddPoint(time, lAdapt.tau(1));
      gr6_3->AddPoint(time, lAdapt.tau(0));

      gra1->AddPoint(time, rAdapt.a_(0));
      gra2->AddPoint(time, rAdapt.a_(1));
      gra3->AddPoint(time, rAdapt.a_(2));
      gra4->AddPoint(time, rAdapt.a_(3));
      gra5->AddPoint(time, rAdapt.a_(4));
      gra6->AddPoint(time, rAdapt.a_(5));
      gra7->AddPoint(time, rAdapt.a_(6));
      gra8->AddPoint(time, rAdapt.a_(7));
      gra9->AddPoint(time, rAdapt.a_(8));
    }

    if(i>2000&&i<10000){
      load << 0,0,-450,0,0,0;
      wk3->setExternalForce(0, load);
    }
    else{ 
      load << 0,0,0,0,0,0;
      wk3->clearExternalForcesAndTorques();
    }

    // if(i>3000&&i<8000){
    //   load << 0,0,-850,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>8000&&i<13000){
    //   load << 0,0,-1050,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else{
    //   load << 0,0,0,0,0,0;
    //   wk3->clearExternalForcesAndTorques();
    // }

    // if(i>3000&&i<6000){
    //   load << 0,0,-50,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>6000&&i<9000){
    //   load << 0,0,-150,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>9000&&i<12000){
    //   load << 0,0,-250,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>12000&&i<15000){
    //   load << 0,0,-350,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>15000&&i<18000){
    //   load << 0,0,-450,0,0,0;
    //   wk3->setExternalForce(0, load);
    // }
    // else if(i>18000&&i<21000){
    //   load << 0,0,-550,0,0,0;
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

    c2->Divide(2,1);
    c2->cd(1);
    mgA1->Add(gra1);
    mgA1->Add(gra2);
    mgA1->Add(gra3);
    mgA1->Add(gra4);
    mgA1->Add(gra5);
    mgA1->Add(gra6);
    mgA1->Draw("AL");
    c2->cd(2);
    mgA2->Add(gra7);
    mgA2->Add(gra8);
    mgA2->Add(gra9);
    mgA2->Draw("AL");
    app.Run();
  }

  server.killServer();
}
