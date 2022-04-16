#ifndef ADAPTIVE_H_
#define ADAPTIVE_H_

#include "kinematics.h"

class WKAdaptive{
    public:
    double thetab, thetah, thetak;
    Vec2 q,dq,ddq;
    Vec2 qd,dqd,ddqd;
    Vec2 qe,dqe,ddqe;
    Vec2 qr,dqr,ddqr;
    Vec2 s;
    Vec2 tau,tau_;
    Vec5 a,a_,da_,ae;
    Mat2 lambda; 
    Mat5 gamma;
    Mat2 Kd; 
    Eigen::Matrix< double , 2 , 5> Y;
    Mat2 H,C,g,H_,C_,g_;

    WKAdaptive();
    void inita(Vec2 _tau);
    void updateq(Vec2 _q, Vec2 _dq, Vec2 _qd, Vec2 _dqd);
    void updatea();
    void updateY();
    void updatetau();
};

#endif