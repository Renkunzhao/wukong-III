#include <iostream>
#include "kinematics.h"

int main()
{
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

    //2022.3.19 01:55
    //雅可比矩阵和正运动学解算完成，在左腿处与学长计算结果一致，在右腿处当hipy有变化时不一致
    //计划用数值方法实现逆运动学算法，请教学长几何法思路
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

    // Eigen::VectorXd pos(28);
    // pos <<  0,0,0,
    //         0,0,0,0,
    //         0,
    //         0,0,0,0,
    //         0,0,0,0,
    //         0,0,0,0,0,0,
    //         0,0,0,0,0,0;
    // Vec3 deltap;
    // deltap << 0,-0.045,0.102;
    // Mat4 T,T_inv;

    // WKLegKinematics wkKin(pos);
    // WKLegKinematicsFoot wkKinL(pos, "LeftLeg");
    // T.setIdentity();
    // T.block(0,0,3,3) = wkKinL.wkLink[6].R;
    // T.block(0,3,3,1) = wkKinL.wkLink[6].p;
    // T_inv = T.inverse();
    // T_inv.block(0,3,3,1) -= deltap;
    // cout << "body: left kinematics " << wkKin.wkLink[12].p.transpose() << endl;
    // cout << "foot: left kinematics " << T_inv.block(0,3,3,1).transpose() << endl;

    // // wkKin.wkLink[7].q = 0.1;
    // wkKin.wkLink[8].q = -0.77;
    // wkKin.wkLink[9].q = -0.92;
    // wkKin.wkLink[10].q = 0.36;
    // wkKin.wkLink[11].q = 0.48;
    // wkKin.wkLink[12].q = pi/2;
    // wkKin.forwardKin(0);
    // cout << "body: left kinematics " << wkKin.wkLink[12].p.transpose() << endl;

    // // wkKinL.wkLink[6].q = -0.1;
    // wkKinL.wkLink[5].q = 0.77;
    // wkKinL.wkLink[4].q = 0.92;
    // wkKinL.wkLink[3].q = -0.36;
    // wkKinL.wkLink[2].q = -0.48;
    // wkKinL.wkLink[1].q = -pi/2;
    // wkKinL.forwardKin(0);
    // T.setIdentity();
    // T.block(0,0,3,3) = wkKinL.wkLink[6].R;
    // T.block(0,3,3,1) = wkKinL.wkLink[6].p;
    // T_inv = T.inverse();
    // T_inv.block(0,3,3,1) -= deltap;
    // cout << "foot: left kinematics " << T_inv.block(0,3,3,1).transpose() << endl;
}
