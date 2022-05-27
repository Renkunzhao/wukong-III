#include <iomanip>
#include <string>
#include <random>
#include <vector>
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

void mean_stdev(vector<double> array){
  double sum = accumulate(begin(array), end(array), 0.0);
	double mean =  sum / array.size(); //均值

	double accum  = 0.0;
	for_each (begin(array), end(array), [&](const double d) {
		accum  += (d-mean)*(d-mean);
	});

	double stdev = sqrt(accum/(array.size()-1)); //方差

  cout << "mean:" << mean << endl;
  cout << "stdev:" << stdev << endl;
}

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  //创建世界、导入wk3
  raisim::World world;
  world.setTimeStep(0.001);   // 仿真环境中的更新间隔
  world.addGround();
  // world.setGravity(Eigen::VectorXd::Zero(3));
  auto g = world.getGravity();
  auto dt = world.getTimeStep();
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
  Vec6 gTorso, fTorso, load, footFix;
  gTorso << (totalMass-mass[15]-mass[21])*g.e(),0,0,0;
  fTorso.setZero();
  load.setZero();
  cout << "gTorso: " << gTorso.transpose() << endl;

  //初始化WKLegKinematics类，根据需要的质心高度计算各个关节角度
  Eigen::VectorXd posInit(posDim), posTarget(posDim), velTarget(dof), feedForwardF(dof);
  posTarget <<  0,0, 0.75, 1,0,0,0,
                0., 
                0,0,0,0, 0,0,0,0,
                0., 0.042578, -0.492781, 0.97759, -0.03982, -0.535,
                0., -0.042567, -0.498771, 0.97759, 0.03964, -0.535;
                // 1.95698e-07,  7.98862e-07,  0.594518,     1.53575,      -1.707e-07,   -1.13027,     
                // -2.90161e-07, -1.18447e-06, 0.594518,     1.53575,      -5.2841e-07,  -1.13027;    //0.,-0.122,-0.6,theta=-1.3
  WKLegKinematics wkKin(posTarget);
  double pitchTheta = 0.1;
  Eigen::AngleAxisd pitchAngleCoM(pitchTheta,Eigen::Vector3d::UnitY());
  Eigen::Quaterniond quaternion(pitchAngleCoM);
  Eigen::AngleAxisd pitchAngle(-pitchTheta,Eigen::Vector3d::UnitY());
  Vec3 endPr, endPl;
  Mat3 endR;
  // endPr << 0.023,-0.122,-0.8;
  // endPl << 0.023,0.122,-0.8;
  endPr << 0.,-0.122,-0.7;
  endPl << 0.,0.122,-0.7;
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
  TMultiGraph *mg2 = new TMultiGraph();
  auto gr2_1 = new TGraph();
  auto gr2_2 = new TGraph();
  TMultiGraph *mg3 = new TMultiGraph();
  auto gr3_1 = new TGraph();
  auto gr3_2 = new TGraph();
  TMultiGraph *mg4 = new TMultiGraph();
  auto gr4_1 = new TGraph();
  auto gr4_2 = new TGraph();
  TMultiGraph *mg5 = new TMultiGraph();
  auto gr5_1 = new TGraph();
  auto gr5_2 = new TGraph();
  TMultiGraph *mg6 = new TMultiGraph();
  auto gr6_1 = new TGraph();
  auto gr6_2 = new TGraph();

  TCanvas* c1 = new TCanvas("c1", "CoM", 0, 0, 1400, 600);
  TMultiGraph *mgCoM = new TMultiGraph();
  auto grCoMz = new TGraph();  
  auto grCoMtheta = new TGraph();  

  //快速站稳
  raisim::Mat<3, 3> bodyRotation;
  Vec3 thetad,theta;
  thetad << 0,pitchTheta,0;
  for (int i=0; i<2000; i++) {
    auto pos = wk3->getGeneralizedCoordinate().e();
    auto vel = wk3->getGeneralizedVelocity().e();
    wk3->getBodyOrientation(0, bodyRotation);
    quaToRpy(bodyRotation, theta);
    wkKin.updateLinkq(pos);

    //根据雅可比矩阵计算关节力矩
    fTorso << 0,0,0*(0.8-pos(2)) + 0*(0-vel(2)),
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

  vector<double> ankleTau, q1, q2, q3;
  for (int i=0; i<10000; i++) {
    auto time = world.getWorldTime();
    auto pos = wk3->getGeneralizedCoordinate().e();
    auto vel = wk3->getGeneralizedVelocity().e();
    auto force = wk3->getGeneralizedForce().e();
    wk3->getBodyOrientation(0, bodyRotation);
    quaToRpy(bodyRotation, theta);
    wkKin.updateLinkq(pos);

    if(i>9000){
      ankleTau.push_back(force(20));
      q1.push_back(pos(21));
      q2.push_back(pos(19));
      q3.push_back(pos(18));
    }

    //根据雅可比矩阵计算关节力矩
    fTorso << 0,0,0*(0.8-pos(2)) + 0*(0-vel(2)),
              0*(thetad[0]-theta[0]) + 0*(0-vel(3)),
              2500*(thetad[1]-theta[1]) + 20*(0-vel(4)),
              0;
    wkKin.RTorque = wkKin.RJac.transpose()*((gTorso-fTorso+load)/2);
    wkKin.LTorque = wkKin.LJac.transpose()*((gTorso-fTorso+load)/2);
    feedForwardF << Eigen::VectorXd::Zero(dof-12), wkKin.RTorque, wkKin.LTorque;
    wk3->setGeneralizedForce(feedForwardF);

    if(i%250==0){

      cout << "worldtime: " << time << endl;
      cout << "thetad: " << thetad.transpose() << endl;
      cout << "posTarget: " << posTarget.tail(12).transpose() << endl;

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
    }

    this_thread::sleep_for(chrono::microseconds(900));
    server.integrateWorldThreadSafe();
  }

  mean_stdev(ankleTau);
  mean_stdev(q1);
  mean_stdev(q2);
  mean_stdev(q3);

  {//ROOT设置Pad等操作
    gr1_2->SetLineColor(2);
    gr2_2->SetLineColor(2);
    gr3_2->SetLineColor(2);
    gr4_2->SetLineColor(2);
    gr5_2->SetLineColor(2);
    gr6_2->SetLineColor(2);

    c1->Divide(2,1);
    c1->cd(1);
    grCoMz->Draw("AL");
    c1->cd(2);
    grCoMtheta->Draw("AL");

    c->Divide(3,2);
    c->cd(1);
    mg1->Add(gr1_1);
    mg1->Add(gr1_2);
    mg1->Draw("AL");
    c->cd(2);
    mg2->Add(gr2_1);
    mg2->Add(gr2_2);
    mg2->Draw("AL");
    c->cd(3);
    mg3->Add(gr3_1);
    mg3->Add(gr3_2);
    mg3->Draw("AL");
    c->cd(4);
    mg4->Add(gr4_1);
    mg4->Add(gr4_2);
    mg4->Draw("AL");
    c->cd(5);
    mg5->Add(gr5_1);
    mg5->Add(gr5_2);
    mg5->Draw("AL");
    c->cd(6);
    mg6->Add(gr6_1);
    mg6->Add(gr6_2);
    mg6->Draw("AL");

    app.Run();
  }

  server.killServer();
}
