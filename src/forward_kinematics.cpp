#include <fstream>
#include <iostream>
#include "kinematics.h"
#include "raisim/RaisimServer.hpp"
#include "raisim/World.hpp"
#include "matplotlibcpp.h"

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
    Eigen::VectorXd pos(28);
    pos <<  0,0,0,
            0,0,0,0,
            0,
            0,0,0,0,
            0,0,0,0,
            0., 0.042578, -0.492781, 0.97759, -0.03982, -0.535,
            0., -0.042567, -0.498771, 0.97759, 0.03964, -0.535;
    WKLegKinematics wkKin(pos);
    // wkKin.wkLink[7].q = 0.1;
    // wkKin.wkLink[8].q = -0.77;
    // wkKin.wkLink[9].q = -0.92;
    // wkKin.wkLink[10].q = 0.35;
    // wkKin.forwardKin(0);
    // for(int i=7;i<=12;i++){
    //     cout << wkKin.wkLink[i].name << endl;
    //     cout << "p: " << wkKin.wkLink[i].p.transpose() << endl;
    //     cout << "R:\n" << wkKin.wkLink[i].R.transpose() << endl;
    // }

    // wkKin.wkLink[1].q = 0.1;
    // wkKin.wkLink[2].q = 1;
    // wkKin.wkLink[3].q = pi/2;
    // wkKin.wkLink[4].q = -0.2;
    // wkKin.forwardKin(0);
    // for(int i=1;i<=6;i++){
    //     cout << wkKin.wkLink[i].name << endl;
    //     cout << "p: " << wkKin.wkLink[i].p.transpose() << endl;
    //     cout << "R:\n" << wkKin.wkLink[i].R.transpose() << endl;
    // }

    // wkKin.Jacbian();
    // cout << "left jacbian\n" << wkKin.LJac << endl;
    // cout << "right jacbian\n" << wkKin.RJac << endl;

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
    gc <<  0,0, 2, 1,0,0,0,
                0., 
                0,0,0,0, 0,0,0,0,
                0., 0.042578, -0.492781, 0.97759, -0.03982, -0.535,
                0., -0.042567, -0.498771, 0.97759, 0.03964, -0.535;
    for(int i=0;i<12;i++){
        gc[i+16] = wkKin.wkLink[i+1].q;
    }
    wk3->setGeneralizedCoordinate(gc);
    wk3->setGeneralizedVelocity(Eigen::VectorXd::Zero(dof));  

    /// launch raisim server
    //  运行前最好先打开unity可视化界面
    raisim::RaisimServer server(&world);
    server.launchServer();
    server.focusOn(wk3);
    std::this_thread::sleep_for(std::chrono::milliseconds(2000));

    double CoMx, CoMz;
    double CoMTheta = 0.;
    Eigen::AngleAxisd pitchAngle(-CoMTheta,Eigen::Vector3d::UnitY());
    Vec6 endVr,endVl;
    Vec3 endPr,endPl;
    Mat3 endR;
    Vec6 dqr,dql;
    vector<double> q1d,q2d,q3d,dq1d,dq2d,dq3d;
    vector<double> pltt,pltzd,pltthetad;
    for (int i=0; i<10000; i++) {
        auto time = world.getWorldTime();
        CoMx = 0;
        CoMz = 0.6 + 0.1*sin(0.5*pi*(time-1));
        endVr << 0,0,-0.1*0.5*pi*cos(0.5*pi*(time-1)),0,0,0;
        endVl << 0,0,-0.1*0.5*pi*cos(0.5*pi*(time-1)),0,0,0;

        if(i%100==0){
            //pos
            endPr << -CoMx,-0.122,-CoMz;
            endPl << -CoMx,0.122,-CoMz;
            endPr = pitchAngle.matrix()*endPr;
            endPl = pitchAngle.matrix()*endPl;
            endR = pitchAngle.matrix();
            wkKin.inverseKin(endPr, endR, "RightLeg");
            wkKin.inverseKin(endPl, endR, "LeftLeg");
            // vel
            dqr = wkKin.RJac.fullPivLu().solve(endVr);
            dql = wkKin.LJac.fullPivLu().solve(endVl);
        }

        if(i%100==0){
            cout << "worldtime: " << time << endl;
            pltt.push_back(time);
            pltzd.push_back(CoMz);
            pltthetad.push_back(CoMTheta);
            q1d.push_back(wkKin.wkLink[6].q);
            q2d.push_back(wkKin.wkLink[4].q);
            q3d.push_back(wkKin.wkLink[3].q);
            dq1d.push_back(dqr[5]);
            dq2d.push_back(dqr[3]);
            dq3d.push_back(dqr[2]);
        }

        for(int i=0;i<12;i++){
            gc[i+16] = wkKin.wkLink[i+1].q;
        }
        wk3->setGeneralizedCoordinate(gc);
        wk3->setGeneralizedVelocity(Eigen::VectorXd::Zero(dof)); 

        this_thread::sleep_for(chrono::microseconds(900));
        server.integrateWorldThreadSafe();
    }

    plt::figure(1);
    plt::plot(pltt, pltzd, "k-",{{"label", "z$_{set}$"}});
    plt::legend();
    plt::xlabel("时间(s)");
    plt::ylabel("高度(m)");
    plt::title("机器人质心高度$z$");

    plt::figure(2);
    plt::plot(pltt, pltthetad, "k-", {{"label", "$\\theta_{set}$"}});
    plt::legend();
    plt::xlabel("时间(s)");
    plt::ylabel("角度(rad)");
    plt::title("机器人俯仰角度$\\theta$");

    plt::figure(3);
    plt::plot(pltt, q3d, "m-", {{"label", "HipY$_{set}$"}});
    plt::legend();
    plt::xlabel("时间(s)");
    plt::ylabel("角度(rad)");
    plt::title("关节角度HipY");

    plt::figure(4);
    plt::plot(pltt, q2d, "g-", {{"label", "KneeY$_{set}$"}});
    plt::legend();
    plt::xlabel("时间(s)");
    plt::ylabel("角度(rad)");
    plt::title("关节角度KneeY");

    plt::figure(5);
    plt::plot(pltt, q1d, "r-", {{"label", "AnkleY$_{set}$"}});
    plt::legend();
    plt::xlabel("时间(s)");
    plt::ylabel("角度(rad)");
    plt::title("关节角度AnkleY");

    plt::figure(6);
    plt::plot(pltt, dq3d, "m-", {{"label", "HipY$_{set}$"}});
    plt::legend();
    plt::xlabel("时间(s)");
    plt::ylabel("角速度(rad/s)");
    plt::title("关节角速度HipY");

    plt::figure(7);
    plt::plot(pltt, dq2d, "g-", {{"label", "KneeY$_{set}$"}});
    plt::legend();
    plt::xlabel("时间(s)");
    plt::ylabel("角速度(rad/s)");
    plt::title("关节角速度KneeY");

    plt::figure(8);
    plt::plot(pltt, dq1d, "r-", {{"label", "AnkleY$_{set}$"}});
    plt::legend();
    plt::xlabel("时间(s)");
    plt::ylabel("角速度(rad/s)");
    plt::title("关节角速度AnkleY");

    plt::show();
}
