#include <fstream>
#include <iostream>
#include "kinematics.h"
#include "raisim/RaisimServer.hpp"
#include "raisim/World.hpp"

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
    Eigen::VectorXd pos(28);
    pos <<  0,0,0,
            0,0,0,0,
            0,
            0,0,0,0,
            0,0,0,0,
            0,0,0,0,0,0,
            0,0,0,0,0,0;
    WKLegKinematics wkKin(pos);
    wkKin.wkLink[7].q = 0.1;
    wkKin.wkLink[8].q = -0.77;
    wkKin.wkLink[9].q = -0.92;
    wkKin.wkLink[10].q = 0.35;
    wkKin.forwardKin(0);
    cout << "r: left kinematics " << wkKin.wkLink[12].p.transpose() << endl;

    Vec3 endPos;
    Vec4 nowLegPos;
    nowLegPos <<  0.1,-0.77,-0.92,0.35;
    Fwd_Kin(endPos, "LeftLeg", nowLegPos);  
    cout << "x: left kinematics " << endPos.transpose() << endl;

    wkKin.wkLink[1].q = 0.1;
    wkKin.wkLink[2].q = 1;
    wkKin.wkLink[3].q = pi/2;
    wkKin.wkLink[4].q = -0.2;
    wkKin.forwardKin(0);
    cout << "r: right kinematics " << wkKin.wkLink[6].p.transpose() << endl;

    nowLegPos <<  0.1,-1,-pi/2,0.2;
    Fwd_Kin(endPos, "RightLeg", nowLegPos);  
    cout << "x: right kinematics " << endPos.transpose() << endl;

    Mat6 RJac,LJac;
    wkKin.Jacbian();
    cout << "r: left jacbian\n" << wkKin.LJac << endl;

    Vec6 nowLegPos1;
    nowLegPos1 <<  0.1,-0.77,-0.92,0.35,0,0;
    Jacbian(LJac, "LeftLeg", nowLegPos1);
    cout << "x: left jacbian\n" << LJac << endl;

    wkKin.wkLink[5].q = 1;
    wkKin.wkLink[6].q = 1.1;
    wkKin.forwardKin(0);
    wkKin.Jacbian();
    cout << "r: right jacbian\n" << wkKin.RJac << endl;

    nowLegPos1 <<  0.1,-1,-pi/2,0.2,-1,-1.1;
    Jacbian(RJac, "RightLeg", nowLegPos1);
    cout << "x: right jacbian\n" << RJac << endl;

    Vec3 endPos1;
    Mat3 endR;
    for(int i=0;i<13;i++)
    {
        wkKin.wkLink[i].q = 0.1;
    }
    endPos1 << 0.5,-0.1332,-0.5055;
    endR = Mat3::Identity();
    wkKin.inverseKin(endPos1, endR, "RightLeg");
    Vec6 q;
    for(int i=0;i<6;i++)
    {
        q(i) = wkKin.wkLink[i+1].q;
    }
    cout << "r: right inverse kinematics" << q.transpose() << endl;

    double ang_hipz = 0;
    Vec3 joint_pos;
    Inv_Kin_Pos(endPos1, "RightLeg", ang_hipz, joint_pos);
    cout << "x: right inverse kinematics " << ang_hipz << " " << joint_pos.transpose() << endl;

    for(int i=0;i<13;i++)
    {
        wkKin.wkLink[i].q = 0.1;
    }
    endPos1 << 0.5,0.1332,-0.5055;
    endR = Mat3::Identity();
    wkKin.inverseKin(endPos1, endR, "LeftLeg");
    for(int i=0;i<6;i++)
    {
        q(i) = wkKin.wkLink[i+7].q;  
    }
    cout << "r: left inverse kinematics " << q.transpose() << endl;

    Inv_Kin_Pos(endPos1, "LeftLeg", ang_hipz, joint_pos);
    cout << "x: left inverse kinematics " << ang_hipz << " " << joint_pos.transpose() << endl;

    auto binaryPath = raisim::Path::setFromArgv(argv[0]);

    //创建世界、导入wk3
    raisim::World world;
    world.setTimeStep(0.001);   // 仿真环境中的更新间隔
    world.addGround();
    world.setGravity(Eigen::VectorXd::Zero(3));
    auto wk3 = world.addArticulatedSystem(binaryPath.getDirectory() + "../rsc/wk3-mpcturn/Wukong3_rsm.urdf");
    wk3->setName("wk3");

    //初始姿态
    int posDim = wk3->getGeneralizedCoordinateDim();
    int dof = wk3->getDOF();
    Eigen::VectorXd gc(posDim);
    gc <<  0,0, 2, 0.995,0,0.1,0,
                0., 
                0,0,0,0, 0,0,0,0,
                0., 0.042578, -0.492781, 0.97759, -0.03982, -0.535,
                0., -0.042567, -0.498771, 0.97759, 0.03964, -0.535;
    wkKin.updateLinkq(gc);
    // wk3->setGeneralizedCoordinate(gc);
    // wk3->setGeneralizedVelocity(Eigen::VectorXd::Zero(dof)); 

    //根据需要的质心高度计算各个关节角度
    Eigen::AngleAxisd pitchAngleCoM(0.2,Eigen::Vector3d::UnitY());
    Eigen::Quaterniond quaternion(pitchAngleCoM);
    Eigen::AngleAxisd pitchAngle(-0.2,Eigen::Vector3d::UnitY());
    Vec3 endPr, endPl;
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
        gc[i+16] = wkKin.wkLink[i+1].q;
    }
    gc[3] = quaternion.coeffs()[3];
    gc[4] = quaternion.coeffs()[0];
    gc[5] = quaternion.coeffs()[1];
    gc[6] = quaternion.coeffs()[2];
    cout << "gc: " << gc.transpose() << endl;
    wk3->setGeneralizedCoordinate(gc);
    wk3->setGeneralizedVelocity(Eigen::VectorXd::Zero(dof)); 

    raisim::Mat<3, 3> bodyRotation;
    Vec3 theta;
    wk3->getBodyOrientation(0, bodyRotation);
    quaToRpy(bodyRotation, theta);
    cout << "theta: " << theta.transpose() << endl;

    cout << "M:\n" << wk3->getMassMatrix() << endl;
    cout << "h:\n" << wk3->getNonlinearities(world.getGravity()) << endl;

    /// launch raisim server
    //  运行前最好先打开unity可视化界面
    raisim::RaisimServer server(&world);
    server.launchServer();
    server.focusOn(wk3);
    std::this_thread::sleep_for(std::chrono::milliseconds(2000));

    while(1){
        this_thread::sleep_for(chrono::microseconds(900));
        server.integrateWorldThreadSafe();
    }
}
