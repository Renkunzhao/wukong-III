#ifndef ADAPTIVE3_H_
#define ADAPTIVE3_H_

#include "kinematics.h"

class WKAdaptive3{
    public:
    double thetab, thetah, thetak;
    Vec3 q,dq,ddq;
    Vec3 qd,dqd,ddqd;
    Vec3 qe,dqe,ddqe;
    Vec3 qr,dqr,ddqr;
    Vec3 s;
    Vec3 tau,tau_;
    Vec9 a,a_,da_,ae;
    Mat3 lambda; 
    Mat9 gamma;
    Mat3 Kd; 
    Eigen::Matrix< double , 3 , 9> Y;
    Mat3 H,C,g,H_,C_,g_;

    WKAdaptive3();
    void inita(Vec2 _tau);
    void updates();
    void updatea();
    void updateY();
    void updateY1();
    void updatetau();
};

#endif