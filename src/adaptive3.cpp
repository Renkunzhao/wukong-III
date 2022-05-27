#include <iostream>
#include "adaptive3.h"

using namespace std;

WKAdaptive3::WKAdaptive3(){
    qd.setZero();
    dqd.setZero();
    ddqd.setZero();
    Kd.setZero();
    Kd.diagonal() << 25, 25, 25;
    lambda.setZero();
    lambda.diagonal() <<  8, 8, 8;
    gamma.setZero();
    gamma.diagonal() << 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.2, 0.2, 0.2;
    // gamma.diagonal() << 0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 3, 3, 3;
    a_ << 3,3,2.5,3,1.9,1.9,91,84,53;
}

void WKAdaptive3::updates(){
    // q,dq,qd,dqd,ddqd
    qe = q - qd;
    dqe = dq - dqd;
    dqr = dqd - lambda*qe;
    ddqr = ddqd - lambda*dqe;
    s = dq - dqr;
    if(s.array().abs().maxCoeff()>0.01){
        updateY();
    }
    // else{
    //     updateY1();
    // }
}

void WKAdaptive3::updateY(){
    Y <<    ddqr[0],    ddqr[0]+ddqr[1],    ddqr[0]+ddqr[1]+ddqr[2],  cos(q[1])*(2*ddqr[0]+ddqr[1])-sin(q[1])*(dqr[0]*dq[1]+dqr[1]*(dq[0]+dq[1])),  
                cos(q[1]+q[2])*(2*ddqr[0]+ddqr[1]+ddqr[2])-sin(q[1]+q[2])*(dqr[0]*(dq[1]+dq[2])+dqr[1]*(dq[0]+dq[1]+dq[2])+dqr[2]*(dq[0]+dq[1]+dq[2])),  
                    cos(q[2])*(2*ddqr[0]+2*ddqr[1]+ddqr[2])-sin(q[2])*(dqr[0]*dq[2]+dqr[1]*dq[2]+dqr[2]*(dq[0]+dq[1]+dq[2])),  
                        -sin(q[0]),   -sin(q[0]+q[1]),  -sin(q[0]+q[1]+q[2]),
            0,          ddqr[0]+ddqr[1],    ddqr[0]+ddqr[1]+ddqr[2],  cos(q[1])*ddqr[0]+sin(q[1])*dqr[0]*dq[0],  
                cos(q[1]+q[2])*ddqr[0]+sin(q[1]+q[2])*dqr[0]*dq[0],
                    cos(q[2])*(2*ddqr[0]+2*ddqr[1]+ddqr[2])-sin(q[2])*(dqr[0]*dq[2]+dqr[1]*dq[2]+dqr[2]*(dq[0]+dq[1]+dq[2])),    
                        0,            -sin(q[0]+q[1]),  -sin(q[0]+q[1]+q[2]),
            0,          0,                  ddqr[0]+ddqr[1]+ddqr[2],  0,  
                cos(q[1]+q[2])*ddqr[0]+sin(q[1]+q[2])*dqr[0]*dq[0],  
                    cos(q[2])*(ddqr[0]+ddqr[1])+sin(q[2])*(dqr[0]+dqr[1])*(dq[0]+dq[1]),    
                        0,            0,                -sin(q[0]+q[1]+q[2]),        
    updatea();
}

void WKAdaptive3::updateY1(){    
    Y <<    0,  0,  0,  0,  0,  0,  -sin(q[0]),   -sin(q[0]+q[1]),  -sin(q[0]+q[1]+q[2]),
            0,  0,  0,  0,  0,  0,  0,            -sin(q[0]+q[1]),  -sin(q[0]+q[1]+q[2]),
            0,  0,  0,  0,  0,  0,  0,            0,                -sin(q[0]+q[1]+q[2]),    
    updatea();
}

void WKAdaptive3::updatea(){
    da_ = - gamma * Y.transpose()*s;
    a_ += da_;
    updatetau();
}

void WKAdaptive3::updatetau(){
    tau_ = Y*a_;
    // tau = tau_ - Kd*s;
    tau = tau_;
}