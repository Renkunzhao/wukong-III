#ifndef KINEMATICS_H_
#define KINEMATICS_H_

#include <iostream>
#include <string>
#include <cmath>
#include "Eigen/Dense"
#include <unsupported/Eigen/MatrixFunctions>        //包含矩阵指数运算函数函数

using namespace std;

#define sign(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

const double pi = 3.1415926535898;

typedef Eigen::Matrix< double , 2 , 1> Vec2;
typedef Eigen::Matrix< double , 3 , 1> Vec3;
typedef Eigen::Matrix< double , 4 , 1> Vec4;
typedef Eigen::Matrix< double , 5 , 1> Vec5;
typedef Eigen::Matrix< double , 6 , 1> Vec6;
typedef Eigen::Matrix< double , 9 , 1> Vec9;
typedef Eigen::Matrix< double , 2 , 2> Mat2;
typedef Eigen::Matrix< double , 3 , 3> Mat3;
typedef Eigen::Matrix< double , 4 , 4> Mat4;
typedef Eigen::Matrix< double , 5 , 5> Mat5;
typedef Eigen::Matrix< double , 6 , 6> Mat6;
typedef Eigen::Matrix< double , 9 , 9> Mat9;

struct uLink
{
  std::string name;
  int sister    =-1;
  int child     =-1;
  int mother    =-1;
  int m;        //质量
  Vec3 a;        //关节轴矢量
  Vec3 b;        //相对位置
  double q;        //关节角    
  Vec3 p;        //世界坐标系中的位置
  Mat3 R;        //世界坐标系中的姿态
};

Mat3 hat(Vec3 w);
Vec3 vee(Mat3 S);
Vec6 getVWerr(uLink target, uLink now);

class WKLegKinematics
{
  private: 
  public:
    uLink wkLink[13];
    Mat6 LJac,RJac;
    Vec6 LTorque,RTorque;
    WKLegKinematics(Eigen::VectorXd gc);
    void initLink();
    void printLinkInfo(int i);
    void updateLinkq(Eigen::VectorXd gc);
    void forwardKin(int j);
    void Jacbian();
    void qRangeHandler();
    void inverseKin(Vec3 endPos, Mat3 endR, const string& LegName);
};

class WKLegKinematicsFoot
{
  private: 
  public:
    uLink wkLink[7];
    Mat6 Jac;
    Vec6 Torque;
    WKLegKinematicsFoot(Eigen::VectorXd gc, const string& LegName);
    void initLink(const string& LegName);
    void printLinkInfo(int i);
    void updateLinkq(Eigen::VectorXd gc, const string& LegName);
    void forwardKin(int j);
    void Jacbian();
    void qRangeHandler();
    void inverseKin(Vec3 endPos, Mat3 endR, const string& LegName);
    void jointTorque(Vec6 gHipy, Vec6 gKneey, Vec6 load);
};

void Fwd_Kin(Vec3& endPos, const string& LegName, Vec4 nowLegPos);
void Inv_Kin_Pos(Vec3 endPos, const string& LegName, double& ang_hipz, Vec3& joint_pos);
void Jacbian(Mat6& Jac, const string& LegName, Vec6 nowLegPos);

#endif