#include <iostream>
#include "adaptive.h"

using namespace std;

WKAdaptive::WKAdaptive(){
    qd      <<  -0.513256,  1.15416;
    dqd     <<  0,          0;
    ddqd    <<  0,          0;
    Kd.setZero();
    // Kd.diagonal() << 1.0, 0.7;
    Kd.diagonal() << 25, 25;
    lambda.setZero();
    lambda.diagonal() <<  16, 10;
    lambda.diagonal() <<  8, 8;
    gamma.setZero();
    gamma.diagonal() << 0.005, 0.001, 0.001, 0.30, 0.08;
    a_.setZero();
    a_ << 0,0,0,100,100;
}

void WKAdaptive::updateq(Vec2 _q, Vec2 _dq, Vec2 _qd, Vec2 _dqd){
    //这里dq是直接用raisim读速度还是用q做差分呢？
    // q[0] = thetab + thetah;
    // q[1] = thetab + thetah + thetak;
    q = _q;
    dq = _dq;
    qd = _qd;
    dqd = _dqd;
    updates();
}

void WKAdaptive::updates(){
    qe = q - qd;
    dqe = dq - dqd;
    dqr = dqd - lambda*qe;
    s = dq - dqr;
    // cout << (Kd*s).array().abs().maxCoeff() << endl;
    updateY();
}

void WKAdaptive::updateY(){
    // Y <<    ddqr[0],    ddq[0]*cos(thetak) + ddq[1]*cos(thetak) + dq[0]*dqr[0]*sin(thetak) - dq[1]*dqr[1]*sin(thetak),  ddqr[1],    sin(q[0]),  sin(q[1]),
    //         0,          ddqr[0]*cos(thetak) +  dq[0]*dqr[0]*sin(thetak),                                                ddqr[1],    0,          sin(q[1]);
    Y <<    0,  0,  0,    -sin(q[0]),   -sin(q[0]+q[1]),
            0,  0,  0,    0,            -sin(q[0]+q[1]);     
    updatea();
}

void WKAdaptive::updatea(){
    da_ = - gamma * Y.transpose()*s;
    a_ += da_;
    updatetau();
}

void WKAdaptive::updatetau(){
    tau_ = Y*a_;
    // tau = tau_ - Kd*s;
    tau = tau_;
}