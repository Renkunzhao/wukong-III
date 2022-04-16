#include <iostream>
#include "kinematics.h"

using namespace std;

WKLegKinematicsFoot::WKLegKinematicsFoot(Eigen::VectorXd gc, const string& LegName)
{
    initLink(LegName);
    printLinkInfo(0);
    updateLinkq(gc, LegName);
    cout << "WKLegKinematicsFoot: 腿部连杆运动学信息初始化完成" << endl;
}

void WKLegKinematicsFoot::initLink(const string& LegName)
{
    wkLink[0].name = "foot";
    wkLink[1].name = "AnkleY";
    wkLink[2].name = "AnkleX";
    wkLink[3].name = "KneeY";
    wkLink[4].name = "HipY";
    wkLink[5].name = "HipX";
    wkLink[6].name = "HipZ";
    
    wkLink[0].child = 1;
    wkLink[1].child = 2;
    wkLink[2].child = 3;
    wkLink[3].child = 4;
    wkLink[4].child = 5;
    wkLink[5].child = 6;
    wkLink[6].child = -1;

    wkLink[0].mother = -1;
    wkLink[1].mother = 0;
    wkLink[2].mother = 1;
    wkLink[3].mother = 2;
    wkLink[4].mother = 3;
    wkLink[5].mother = 4;
    wkLink[6].mother = 5;

    wkLink[0].a << 0,0,0;
    wkLink[1].a << 0,1,0;
    wkLink[2].a << 1,0,0;
    wkLink[3].a << 0,1,0;
    wkLink[4].a << 0,1,0;
    wkLink[5].a << 1,0,0;
    wkLink[6].a << 0,0,1;

    if(LegName == "RightLeg"){
        wkLink[0].b << 0,0,0;
        // wkLink[1].b << 0,0,0.039;
        wkLink[1].b << 0,0,0;
        wkLink[2].b << 0,0,0;
        wkLink[3].b << -0.023,0,0.35;
        wkLink[4].b << 0,0,0.35;
        wkLink[5].b << 0,0.03,0;
        wkLink[6].b << 0,0.047,0.0985;
    }
    else if(LegName == "LeftLeg"){
        wkLink[0].b << 0,0,0;
        // wkLink[1].b << 0,0,0.039;
        wkLink[1].b << 0,0,0;
        wkLink[2].b << 0,0,0;
        wkLink[3].b << -0.023,0,0.35;
        wkLink[4].b << 0,0,0.35;
        wkLink[5].b << 0,-0.03,0;
        wkLink[6].b << 0,-0.047,0.0985;
    }
    else{
        cout << "LegName error!!!" << endl;
    }

    wkLink[0].p << 0,0,0;
    wkLink[0].R << 1,0,0,
                    0,1,0,
                    0,0,1;
}

void WKLegKinematicsFoot::printLinkInfo(int i)
{
    if(i!=-1){
        cout << wkLink[i].name <<endl;
        printLinkInfo(wkLink[i].sister);
        printLinkInfo(wkLink[i].child);
    }
}

void WKLegKinematicsFoot::updateLinkq(Eigen::VectorXd gc, const string& LegName)
{
    if(LegName == "RightLeg"){
        for(int i=1;i<=6;i++)
        {
            wkLink[i].q = -gc(22-i);
        }
    }
    else if(LegName == "LeftLeg"){
        for(int i=1;i<=6;i++)
        {
            wkLink[i].q = -gc(28-i);
        }
    }
    else{
        cout << "LegName error!!!" << endl;
    }
    cout << "WKLegKinematicsFoot: 关节角更新完成" << endl;
    forwardKin(0);
    cout << "WKLegKinematicsFoot: 正运动学计算完成" << endl;
    Jacbian();
    cout << "WKLegKinematicsFoot: 雅可比矩阵计算完成" << endl;
}

void WKLegKinematicsFoot::forwardKin(int j)
{
    if(j==-1) return;
    if(j!=0)
    {
        int i = wkLink[j].mother;
        wkLink[j].p = wkLink[i].p + wkLink[i].R * wkLink[j].b;
        wkLink[j].R = wkLink[i].R * (hat(wkLink[j].a)*wkLink[j].q).exp(); 
    }
    forwardKin(wkLink[j].sister);
    forwardKin(wkLink[j].child);
}

void WKLegKinematicsFoot::Jacbian()
{
    forwardKin(0);
    Vec3 a,end;
    int j;
    int idx[6];

    //计算右腿雅可比矩阵
    for(int n=0;n<6;n++){
        idx[n] = n+1;
    }
    end = wkLink[idx[5]].p;
    for(int n=0;n<6;n++){
        j = idx[n];
        a = wkLink[j].R*wkLink[j].a;
        Jac.col(n) <<   a.cross(end-wkLink[j].p), 
                        a;
    }
}

void WKLegKinematicsFoot::qRangeHandler()                    //不同关节有不同关节角范围，这里暂时全部设为+-pi
{
    for(int i=0;i<13;i++)
    {
        wkLink[i].q = fmod(wkLink[i].q, pi);
    }
}

void WKLegKinematicsFoot::inverseKin(Vec3 endPos, Mat3 endR, const string& LegName)
{
    // uLink targetLink;
    // targetLink.p = endPos;
    // targetLink.R = endR;
    // double lambda = 0.5;
    // forwardKin(0);
    // int idx[6];
    // if(LegName == "RightLeg"){
    //     for(int n=0;n<6;n++){
    //         idx[n] = n+1;
    //     }
    // }
    // else if(LegName == "LeftLeg"){
    //     for(int n=0;n<6;n++){
    //         idx[n] = n+7;
    //     }
    // }
    // else{
    //     cout << "LegName error!!!" << endl;
    // }
    // Mat6 Jac;
    // Vec6 err,dq;
    // int n,j;
    // for(n=0;n<1000;n++){
    //     Jacbian(Jac, LegName);
    //     err = getVWerr(targetLink, wkLink[idx[5]]);
    //     if(err.norm()<1e-6){
    //         qRangeHandler();
    //         cout << "InverseKinematics complete" << endl;
    //         cout << "iteration number: " << n << endl;
    //         cout << "error: " << err.norm() << endl;
    //         return;
    //     } 
    //     dq = lambda * (Jac.fullPivLu().solve(err));     //为防止奇异姿态下雅可比不可逆，采用LU分解，但是如果初始姿态为奇异姿态，貌似无法向真值靠近
    //     for(int nn=0;nn<6;nn++){
    //         j = idx[nn];
    //         wkLink[j].q = wkLink[j].q + dq(nn); 
    //     }
    //     forwardKin(0);
    // }
    // cout << "InverseKinematics error!!!" << endl;
    // cout << "iteration number: " << n << endl;
    // cout << "error: " << err.norm() << endl;
}
  
void WKLegKinematicsFoot::jointTorque(Vec6 gHipy, Vec6 gKneey, Vec6 load)
{
    auto rTHipy = Jac.transpose()*(load/2);
    auto rTKneey = Jac.transpose()*((gKneey+load)/2);
    Torque <<  0, 0, rTHipy[2], rTKneey[3], 0, rTHipy[5];
    // RTorque <<  0, 0, 0, rTKneey[3], 0, 0;

    cout << "rTHipy: " << rTHipy.transpose() << endl;
}

