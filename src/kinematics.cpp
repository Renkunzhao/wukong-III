#include <iostream>
#include "kinematics.h"

using namespace std;

Mat3 hat(Vec3 w)
{
    Mat3 S;
    S<< 0,      -w[2],  w[1],
        w[2],   0,      -w[0],
        -w[1],  w[0],   0;     
    return S;
}

Vec3 vee(Mat3 S)
{
    Vec3 w;
    w<< S(2,1), S(0,2), S(1,0);
    return w;
}

Vec6 getVWerr(uLink target, uLink now)
{
    Mat3 dR;
    Vec3 dw;
    Vec6 err;
    dR = now.R.transpose()*target.R;
    dw = vee(dR.log());
    err <<  (target.p-now.p),      //位置偏差
            dw;                    //姿态误差，暂时设为零，之后可以用来计算用于保证脚板水平的踝关节关节角    
    return err;
}

WKLegKinematics::WKLegKinematics(Eigen::VectorXd gc)
{
    initLink();
    printLinkInfo(0);
    updateLinkq(gc);
    cout << "WKLegKinematics: 腿部连杆运动学信息初始化完成" << endl;
}

void WKLegKinematics::initLink()
{
    wkLink[0].name = "body";
    wkLink[1].name = "rHipZ";
    wkLink[2].name = "rHipX";
    wkLink[3].name = "rHipY";
    wkLink[4].name = "rKneeY";
    wkLink[5].name = "rAnkleX";
    wkLink[6].name = "rAnkleY";
    wkLink[7].name = "lHipZ";
    wkLink[8].name = "lHipX";
    wkLink[9].name = "lHipY";
    wkLink[10].name = "lKneeY";
    wkLink[11].name = "lAnkleX";
    wkLink[12].name = "lAnkleY";

    wkLink[1].sister = 7;

    wkLink[0].child = 1;
    wkLink[1].child = 2;
    wkLink[2].child = 3;
    wkLink[3].child = 4;
    wkLink[4].child = 5;
    wkLink[5].child = 6;
    wkLink[6].child = -1;
    wkLink[7].child = 8;
    wkLink[8].child = 9;
    wkLink[9].child = 10;
    wkLink[10].child = 11;
    wkLink[11].child = 12;
    wkLink[12].child = -1;   

    wkLink[0].mother = -1;
    wkLink[1].mother = 0;
    wkLink[2].mother = 1;
    wkLink[3].mother = 2;
    wkLink[4].mother = 3;
    wkLink[5].mother = 4;
    wkLink[6].mother = 5;
    wkLink[7].mother = 0;
    wkLink[8].mother = 7;
    wkLink[9].mother = 8;
    wkLink[10].mother = 9;
    wkLink[11].mother = 10;
    wkLink[12].mother = 11; 

    wkLink[0].a << 0,0,0;
    wkLink[1].a << 0,0,1;
    wkLink[2].a << 1,0,0;
    wkLink[3].a << 0,1,0;
    wkLink[4].a << 0,1,0;
    wkLink[5].a << 1,0,0;
    wkLink[6].a << 0,1,0;
    wkLink[7].a << 0,0,1;
    wkLink[8].a << 1,0,0;
    wkLink[9].a << 0,1,0;
    wkLink[10].a << 0,1,0;
    wkLink[11].a << 1,0,0;
    wkLink[12].a << 0,1,0;

    // wkLink[0].b << 0,0,0;
    // wkLink[1].b << 0,-0.045,-0.0785;
    // wkLink[2].b << 0,-0.047,-0.107;
    // wkLink[3].b << 0,-0.0412,0;
    // wkLink[4].b << 0,0,-0.35;
    // wkLink[5].b << 0,0,-0.35;
    // wkLink[6].b << 0,0,0;
    // wkLink[7].b << 0,0.045,-0.0785;
    // wkLink[8].b << 0,0.047,-0.107;
    // wkLink[9].b << 0,0.0412,0;
    // wkLink[10].b << 0,0,-0.35;
    // wkLink[11].b << 0,0,-0.35;
    // wkLink[12].b << 0,0,0;

    wkLink[0].b << 0,0,0;
    wkLink[1].b << 0,-0.045,-0.102;
    wkLink[2].b << 0,-0.047,-0.0985;
    wkLink[3].b << 0,-0.03,0;
    wkLink[4].b << 0,0,-0.35;
    wkLink[5].b << 0.023,0,-0.35;
    wkLink[6].b << 0,0,0;
    wkLink[7].b << 0,0.045,-0.102;
    wkLink[8].b << 0,0.047,-0.0985;
    wkLink[9].b << 0,0.03,0;
    wkLink[10].b << 0,0,-0.35;
    wkLink[11].b << 0.023,0,-0.35;
    wkLink[12].b << 0,0,0;

    wkLink[0].p << 0,0,0;
    wkLink[0].R <<  1,0,0,
                    0,1,0,
                    0,0,1;
}

void WKLegKinematics::printLinkInfo(int i)
{
    if(i!=-1){
        cout << wkLink[i].name <<endl;
        printLinkInfo(wkLink[i].sister);
        printLinkInfo(wkLink[i].child);
    }
}

void WKLegKinematics::updateLinkq(Eigen::VectorXd gc)
{   
    for(int i=1;i<=12;i++)
    {
        wkLink[i].q = gc(i+15);
    }
    // cout << "WKLegKinematics: 关节角更新完成" << endl;
    forwardKin(0);
    // cout << "WKLegKinematics: 正运动学计算完成" << endl;
    Jacbian();
    // cout << "WKLegKinematics: 雅可比矩阵计算完成" << endl;
}

void WKLegKinematics::forwardKin(int j)
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

void WKLegKinematics::Jacbian()
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
        RJac.col(n) <<  a.cross(end-wkLink[j].p), 
                        a;
    }

    //计算左腿雅可比矩阵
    for(int n=0;n<6;n++){
        idx[n] = n+7;
    }
    end = wkLink[idx[5]].p;
    for(int n=0;n<6;n++){
        j = idx[n];
        a = wkLink[j].R*wkLink[j].a;
        LJac.col(n) <<  a.cross(end-wkLink[j].p), 
                        a;
    }
}

void WKLegKinematics::qRangeHandler()                    //不同关节有不同关节角范围，这里暂时全部设为+-pi
{
    for(int i=0;i<13;i++)
    {
        wkLink[i].q = fmod(wkLink[i].q, pi);
    }
}

void WKLegKinematics::inverseKin(Vec3 endPos, Mat3 endR, const string& LegName)
{
    uLink targetLink;
    targetLink.p = endPos;
    targetLink.R = endR;
    double lambda = 0.5;
    forwardKin(0);
    int idx[6];
    if(LegName == "RightLeg"){
        for(int n=0;n<6;n++){
            idx[n] = n+1;
        }
    }
    else if(LegName == "LeftLeg"){
        for(int n=0;n<6;n++){
            idx[n] = n+7;
        }
    }
    else{
        cout << "LegName error!!!" << endl;
    }
    Mat6 Jac;
    Vec6 err,dq;
    int n,j;
    for(n=0;n<1000;n++){
        Jacbian();
        if(LegName == "RightLeg"){
            Jac = RJac;
        }
        else if(LegName == "LeftLeg"){
            Jac = LJac;
        }
        else{
            cout << "LegName error!!!" << endl;
        }
        err = getVWerr(targetLink, wkLink[idx[5]]);
        if(err.norm()<1e-6){
            qRangeHandler();
            cout << "InverseKinematics complete" << endl;
            cout << "iteration number: " << n << endl;
            cout << "error: " << err.norm() << endl;
            return;
        } 
        dq = lambda * (Jac.fullPivLu().solve(err));     //为防止奇异姿态下雅可比不可逆，采用LU分解，但是如果初始姿态为奇异姿态，貌似无法向真值靠近
        for(int nn=0;nn<6;nn++){
            j = idx[nn];
            wkLink[j].q = wkLink[j].q + dq(nn); 
        }
        forwardKin(0);
    }
    cout << "InverseKinematics error!!!" << endl;
    cout << "iteration number: " << n << endl;
    cout << "error: " << err.norm() << endl;
}

void Fwd_Kin(Vec3& endPos, const string& LegName, Vec4 nowLegPos)
{
  static double posx,posy,posz;
  double t1,t2,t3,t4;
  double c1,c2,c3,c4,s1,s2,s3,s4;
  double delta1,delta2,delta3,delta4,delta5,delta6,delta7;
  Vec3 Hip_Vec, Hip_X_To_Z, Hip_Y_To_X, Knee_To_Hip_Y, EndPoint_To_Knee;

//   t1 = Amp_Hip_X->Angle_Measure();
//   t2 = Amp_Hip_Y->Angle_Measure();
//   t3 = Amp_Hip_Z->Angle_Measure();
//   t4 = Amp_Knee->Angle_Measure();
  t1 = nowLegPos[1];
  t2 = nowLegPos[2];
  t3 = nowLegPos[0];
  t4 = nowLegPos[3];

  c1 = cos(t1);
  s1 = sin(t1);
  c2 = cos(t2);
  s2 = sin(t2);
  c3 = cos(t3);
  s3 = sin(t3);
  c4 = cos(t4);
  s4 = sin(t4);

  double x,y,z,y1,z1,x2,y2,l1,l2;
  
  if(LegName == "LeftLeg"){
    //由机器人各连杆长度表示的，笔直站立状态下各局部坐标系原点平移向量
    Hip_Vec << 0,0.045,-0.0785;
    Hip_X_To_Z << 0,0.047,-0.107;
    Hip_Y_To_X << 0,0.0412,0;
    Knee_To_Hip_Y <<0,0,-0.35;
    EndPoint_To_Knee <<0,0,-0.35;

    x = Hip_Vec[0];
    y = Hip_Vec[1];
    z = Hip_Vec[2];
    y1 = Hip_X_To_Z[1];
    z1 = Hip_X_To_Z[2];
    x2 = Hip_Y_To_X[0];
    y2 = Hip_Y_To_X[1];
    l1 = abs(Knee_To_Hip_Y[2]);
    l2 = abs(EndPoint_To_Knee[2]);
    delta1 = s2*s3-c2*c3*s1;
    delta2 = c3*s2+c2*s1*s3;

    posx = x+x2*c3-y1*s3-l2*(c4*delta2
        +s4*(c2*c3-s1*s2*s3))-l1*delta2-y2*c1*s3;

    posy = y+y1*c3+x2*s3-l2*(c4*delta1
        +s4*(c2*s3+c3*s1*s2))-l1*delta1+y2*c1*c3;

    posz = z+z1-l2*(c1*c2*c4-c1*s2*s4)+y2*s1
         -l1*c1*c2;
  }
  else if(LegName == "RightLeg"){
    Hip_Vec << 0,-0.045,-0.0785;
    Hip_X_To_Z << 0,-0.047,-0.107;
    Hip_Y_To_X << 0,-0.0412,0;
    Knee_To_Hip_Y <<0,0,-0.35;
    EndPoint_To_Knee <<0,0,-0.35;
    x = Hip_Vec[0];
    y = Hip_Vec[1];
    z = Hip_Vec[2];
    y1 = Hip_X_To_Z[1];
    z1 = Hip_X_To_Z[2];
    x2 = Hip_Y_To_X[0];
    y2 = Hip_Y_To_X[1];
    l1 = abs(Knee_To_Hip_Y[2]);
    l2 = abs(EndPoint_To_Knee[2]);
    delta1 = s2*s3-c2*c3*s1;
    delta2 = c3*s2+c2*s1*s3;

    posx = x+x2*c3+y1*s3-l2*(c4*delta2
        +s4*(c2*c3-s1*s2*s3))-l1*delta2+y2*c1*s3;

    posy = y+y1*c3-x2*s3+l2*(c4*delta1
        +s4*(c2*s3+c3*s1*s2))+l1*delta1+y2*c1*c3;

    posz = z+z1-l2*(c1*c2*c4-c1*s2*s4)-y2*s1
         -l1*c1*c2;
    
  }
	endPos << posx, posy, posz;
}

void Inv_Kin_Pos(Vec3 endPos, const string& LegName, double& ang_hipz, Vec3& joint_pos)
{
  Vec3 Hip_Vec, Hip_X_To_Z, Hip_Y_To_X, Knee_To_Hip_Y, EndPoint_To_Knee;
  double x,y,z;
  double l0,l1,l2;
  static double t1,t2,t3;
  const double pi = 3.141593;
  Vec3 Hip_X_To_Robot;
  //double theta_z = Amp_Hip_Z->Angle_Measure();
  double theta_z = ang_hipz;
  Mat3 HipZ_Mat;
//   double x_Hip_y_to_x = abs(Hip_Y_To_X[0]);
//   l0 = abs(Hip_Y_To_X[1]);
//   l1 = abs(Knee_To_Hip_Y[2]);
//   l2 = abs(EndPoint_To_Knee[2]);
  double x_Hip_y_to_x = 0;
  l0 = 0.0412;
  l1 = 0.35;
  l2 = 0.35;
  Vec3 Leg_Vec;

  if(LegName == "LeftLeg"){
    Hip_Vec << 0,0.045,-0.0785;
    Hip_X_To_Z << 0,0.047,-0.107;
    Hip_Y_To_X << 0,0.0412,0;
    Knee_To_Hip_Y <<0,0,-0.35;
    EndPoint_To_Knee <<0,0,-0.35;

    HipZ_Mat<<cos(theta_z),-sin(theta_z),0,
                   sin(theta_z),cos(theta_z),0,
                   0,0,1;
    Hip_X_To_Robot = Hip_Vec + HipZ_Mat*Hip_X_To_Z  ;

    Leg_Vec = HipZ_Mat.transpose()*(endPos-Hip_X_To_Robot);

    x=Leg_Vec[0];
    y=-Leg_Vec[1];
    z=Leg_Vec[2];
    double theta_t = acos(-y/sqrt(y*y+z*z));
    t1 = sign(z)*(-acos(l0/sqrt(y*y + z*z)) + theta_t);
    double la = sqrt(x*x+(-z+l0*sin(t1))*(-z+l0*sin(t1))+(-y-l0*cos(t1))*(-y-l0*cos(t1)));
    t3 = pi - acos((l1*l1+l2*l2-la*la)/(2*l1*l2));
    if(z>l0*cos(t1)) t2 = -asin(l2/la*sin(t3))-(acos(x/la)+pi/2);
    else t2 = -asin(l2/la*sin(t3))-asin(x/la);
  }
  else if(LegName == "RightLeg"){
    Hip_Vec << 0,-0.045,-0.0785;
    Hip_X_To_Z << 0,-0.047,-0.107;
    Hip_Y_To_X << 0,-0.0412,0;
    Knee_To_Hip_Y <<0,0,-0.35;
    EndPoint_To_Knee <<0,0,-0.35;

    HipZ_Mat<<cos(-theta_z),-sin(-theta_z),0,
                   sin(-theta_z),cos(-theta_z),0,
                   0,0,1;
    Hip_X_To_Robot = Hip_Vec + HipZ_Mat*Hip_X_To_Z  ;

    Leg_Vec = HipZ_Mat.transpose()*(endPos-Hip_X_To_Robot);

    x=Leg_Vec[0];
    y=Leg_Vec[1];
    z=Leg_Vec[2];
    double theta_t = acos(-y/sqrt(y*y+z*z));

    t1 = sign(z)*(-acos(l0/sqrt(y*y + z*z)) + theta_t);
    double la = sqrt(x*x+(-z+l0*sin(t1))*(-z+l0*sin(t1))+(-y-l0*cos(t1))*(-y-l0*cos(t1)));
    t3 = pi - acos((l1*l1+l2*l2-la*la)/(2*l1*l2));
    if(z>l0*cos(t1)) t2 = -asin(l2/la*sin(t3))-(acos(x/la)+pi/2);
    else t2 = -asin(l2/la*sin(t3))-asin(x/la);
  }
  else{
    std::cout<<"Leg Name Error!!! Inv Kin Error"<<endl;
    //errorFlag = 1;
  }
  if(isnan(t1)||isnan(t2)||isnan(t3)){
    cout<<"Leg Inv Error!!!"<<t1<<" "<<t2<<" "<<t3<<endl;
    //errorFlag = 1;
  }

  joint_pos<<t1, t2,t3;
}

void Jacbian(Mat6& Jac, const string& LegName, Vec6 nowLegPos)
{
    double c1,c2,c3,c4,c5,s1,s2,s3,s4,s5;
    double l0,l1,l2,l3,l4,l5;
    Vec3 Hip_Vec, Hip_X_To_Z, Hip_Y_To_X, Knee_To_Hip_Y, EndPoint_To_Knee;
    Mat6 Jacbian;

    if(LegName == "LeftLeg")
    {
        Hip_Vec << 0,0.045,-0.0785;
        Hip_X_To_Z << 0,0.047,-0.107;
        Hip_Y_To_X << 0,0.0412,0;
        Knee_To_Hip_Y <<0,0,-0.35;
        EndPoint_To_Knee <<0,0,-0.35;

        l0 = abs(Hip_Vec[2] + Hip_X_To_Z[2]);
        l1 = Hip_Vec[1];
        l2 = Hip_X_To_Z[1];
        l3 = Hip_Y_To_X[1];
        l4 = abs(Knee_To_Hip_Y[2]);
        l5 = abs(EndPoint_To_Knee[2]);
        s1 = sin(nowLegPos[0]); s2 = sin(nowLegPos[1]); s3 = sin(nowLegPos[2]);
        s4 = sin(nowLegPos[3]); s5 = sin(nowLegPos[4]); 
        c1 = cos(nowLegPos[0]); c2 = cos(nowLegPos[1]); c3 = cos(nowLegPos[2]);
        c4 = cos(nowLegPos[3]); c5 = cos(nowLegPos[4]); 

        Jacbian << l4 * (s1 * s3 - c1 * c3 * s2) + l5 * (c4 * (s1 * s3 - c1 * c3 * s2) + s4 * (c3 * s1 + c1 * s2 * s3)) - l2 * c1 - l3 * c1 * c2, l3 * s1 * s2 - l5 * (c2 * c3 * c4 * s1 - c2 * s1 * s3 * s4) - l4 * c2 * c3 * s1, -l4 * (c1 * c3 - s1 * s2 * s3) - l5 * (c4 * (c1 * c3 - s1 * s2 * s3) - s4 * (c1 * s3 + c3 * s1 * s2)), -l5 * (c4 * (c1 * c3 - s1 * s2 * s3) - s4 * (c1 * s3 + c3 * s1 * s2)), 0, 0,
            -l4 * (c1 * s3 + c3 * s1 * s2) - l5 * (c4 * (c1 * s3 + c3 * s1 * s2) + s4 * (c1 * c3 - s1 * s2 * s3)) - l2 * s1 - l3 * c2 * s1, l5 * (c1 * c2 * c3 * c4 - c1 * c2 * s3 * s4) - l3 * c1 * s2 + l4 * c1 * c2 * c3, -l4 * (c3 * s1 + c1 * s2 * s3) - l5 * (c4 * (c3 * s1 + c1 * s2 * s3) - s4 * (s1 * s3 - c1 * c3 * s2)), -l5 * (c4 * (c3 * s1 + c1 * s2 * s3) - s4 * (s1 * s3 - c1 * c3 * s2)), 0, 0,
            0, l3 * c2 - l5 * (s2 * s3 * s4 - c3 * c4 * s2) + l4 * c3 * s2, l5 * (c2 * c3 * s4 + c2 * c4 * s3) + l4 * c2 * s3, l5 * (c2 * c3 * s4 + c2 * c4 * s3), 0, 0,
            0, c1, -c2 * s1, -c2 * s1, c4 * (c1 * c3 - s1 * s2 * s3) - s4 * (c1 * s3 + c3 * s1 * s2), s5 * (c4 * (c1 * s3 + c3 * s1 * s2) + s4 * (c1 * c3 - s1 * s2 * s3)) - c2 * c5 * s1,
            0, s1, c1 * c2, c1 * c2, c4 * (c3 * s1 + c1 * s2 * s3) - s4 * (s1 * s3 - c1 * c3 * s2), s5 * (c4 * (s1 * s3 - c1 * c3 * s2) + s4 * (c3 * s1 + c1 * s2 * s3)) + c1 * c2 * c5,
            1, 0, s2, s2, -c2 * c3 * s4 - c2 * c4 * s3, c5 * s2 + s5 * (c2 * c3 * c4 - c2 * s3 * s4); 
    }
    else if(LegName == "RightLeg")
    {
        Hip_Vec << 0,-0.045,-0.0785;
        Hip_X_To_Z << 0,-0.047,-0.107;
        Hip_Y_To_X << 0,-0.0412,0;
        Knee_To_Hip_Y <<0,0,-0.35;
        EndPoint_To_Knee <<0,0,-0.35;

        l0 = abs(Hip_Vec[2] + Hip_X_To_Z[2]);
        l1 = Hip_Vec[1];
        l2 = Hip_X_To_Z[1];
        l3 = Hip_Y_To_X[1];
        l4 = abs(Knee_To_Hip_Y[2]);
        l5 = abs(EndPoint_To_Knee[2]);
        s1 = sin(-nowLegPos[0]); s2 = sin(-nowLegPos[1]); s3 = sin(nowLegPos[2]);
        s4 = sin(nowLegPos[3]); s5 = sin(-nowLegPos[4]); 
        c1 = cos(-nowLegPos[0]); c2 = cos(-nowLegPos[1]); c3 = cos(nowLegPos[2]);
        c4 = cos(nowLegPos[3]); c5 = cos(-nowLegPos[4]);

        Jacbian << l4 * (s1 * s3 - c1 * c3 * s2) + l5 * (c4 * (s1 * s3 - c1 * c3 * s2) + s4 * (c3 * s1 + c1 * s2 * s3)) - l2 * c1 - l3 * c1 * c2, l3 * s1 * s2 - l5 * (c2 * c3 * c4 * s1 - c2 * s1 * s3 * s4) - l4 * c2 * c3 * s1, -l4 * (c1 * c3 - s1 * s2 * s3) - l5 * (c4 * (c1 * c3 - s1 * s2 * s3) - s4 * (c1 * s3 + c3 * s1 * s2)), -l5 * (c4 * (c1 * c3 - s1 * s2 * s3) - s4 * (c1 * s3 + c3 * s1 * s2)), 0, 0,
            -l4 * (c1 * s3 + c3 * s1 * s2) - l5 * (c4 * (c1 * s3 + c3 * s1 * s2) + s4 * (c1 * c3 - s1 * s2 * s3)) - l2 * s1 - l3 * c2 * s1, l5 * (c1 * c2 * c3 * c4 - c1 * c2 * s3 * s4) - l3 * c1 * s2 + l4 * c1 * c2 * c3, -l4 * (c3 * s1 + c1 * s2 * s3) - l5 * (c4 * (c3 * s1 + c1 * s2 * s3) - s4 * (s1 * s3 - c1 * c3 * s2)), -l5 * (c4 * (c3 * s1 + c1 * s2 * s3) - s4 * (s1 * s3 - c1 * c3 * s2)), 0, 0,
            0, l3 * c2 - l5 * (s2 * s3 * s4 - c3 * c4 * s2) + l4 * c3 * s2, l5 * (c2 * c3 * s4 + c2 * c4 * s3) + l4 * c2 * s3, l5 * (c2 * c3 * s4 + c2 * c4 * s3), 0, 0,
            0, c1, -c2 * s1, -c2 * s1, c4 * (c1 * c3 - s1 * s2 * s3) - s4 * (c1 * s3 + c3 * s1 * s2), s5 * (c4 * (c1 * s3 + c3 * s1 * s2) + s4 * (c1 * c3 - s1 * s2 * s3)) - c2 * c5 * s1,
            0, s1, c1 * c2, c1 * c2, c4 * (c3 * s1 + c1 * s2 * s3) - s4 * (s1 * s3 - c1 * c3 * s2), s5 * (c4 * (s1 * s3 - c1 * c3 * s2) + s4 * (c3 * s1 + c1 * s2 * s3)) + c1 * c2 * c5,
            1, 0, s2, s2, -c2 * c3 * s4 - c2 * c4 * s3, c5 * s2 + s5 * (c2 * c3 * c4 - c2 * s3 * s4); 
    }
    Jac << Jacbian;

}
