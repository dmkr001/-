#ifndef DISTORTION_CAMERA_H
#define DISTORTION_CAMERA_H

#include <generic_camera.h>
#include <visp/vpFeaturePoint3D.h>
#include<iostream>
namespace covis
{

// distortion camera class, from ICRA 2016 paper
// xi = (px, py, u0, v0, alpha, beta)
class DistortionCamera : public GenericCamera
{
public:

    vpFeaturePoint3D p_;
    vpMatrix dPdX_;
    int camera_model;

    DistortionCamera(const double &_px, const double &_py, const double &_u0, const double &_v0, const double &_alpha)
    {
        //UCM
        dPdX_.resize(2,3);
        xi_.resize(5);
        xi_[0] = _px;
        xi_[1] = _py;
        xi_[2] = _u0;
        xi_[3] = _v0;
        xi_[4] = _alpha;
        K.resize(5);
        K[0] = K[1] = K[2] = K[3] = K[4] = 0;
        camera_model=0;
    }

    DistortionCamera(const double &_px, const double &_py, const double &_u0, const double &_v0, const double &_alpha, const double &_beta)
    {
        //EUCM
        dPdX_.resize(2,3);
        xi_.resize(6);
        xi_[0] = _px;
        xi_[1] = _py;
        xi_[2] = _u0;
        xi_[3] = _v0;
        xi_[4] = _alpha;
        xi_[5] = _beta;
        K.resize(5);
        K[0] = K[1] = K[2] = K[3] = K[4] = 0;
        camera_model=1;
    }

    DistortionCamera(const double &_px, const double &_py, const double &_u0, const double &_v0, const double &_alpha, const double &_xi, const int _cm)
    {
        //EUCM
        dPdX_.resize(2,3);
        xi_.resize(6);
        xi_[0] = _px;
        xi_[1] = _py;
        xi_[2] = _u0;
        xi_[3] = _v0;
        xi_[4] = _alpha;
        xi_[5] = _xi;
        K.resize(5);
        K[0] = K[1] = K[2] = K[3] = K[4] = 0;
        camera_model=_cm;
    }

    // compute pixel coordinates of a 3D point
    // we assume the point is already in the camera frame
    void project(const vpPoint &_P, double &_u, double &_v)
    {
        switch(camera_model) {
          case 0:
            // UCM
           {

            const double rho = sqrt(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y()+_P.get_Z()*_P.get_Z());
            double x =(1-xi_[4])*(1-xi_[4])*_P.get_X()/(_P.get_Z()*(1-xi_[4])+rho*xi_[4]);
            double y =(1-xi_[4])*(1-xi_[4])*_P.get_Y()/(_P.get_Z()*(1-xi_[4])+rho*xi_[4]);
            const double d = x*x+y*y;
            const double gamma1=xi_[0]/(1-xi_[4]);
            const double gamma2=xi_[1]/(1-xi_[4]);

            const double L = 1+K[0]*d+K[1]*d*d+K[4]*d*d*d;
            x = x*L+2*K[2]*x*y+K[3]*(d+2*x*x);
            y = y*L+2*K[2]*(d+2*y*y)+K[3]*x*y;


//            const double nu_inv = 1./(xi_[4]*rho + (1-xi_[4])*_P.get_Z());
            _u = gamma1*x + xi_[2];    // u = px.X/nu + u0
            _v = gamma2*y + xi_[3];    // v = py.Y/nu + v0
//            const double rho = sqrt(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y()+_P.get_Z()*_P.get_Z());
//            const double nu_inv = 1./(xi_[4]*rho + (1-xi_[4])*_P.get_Z());
//            _u = xi_[0]*_P.get_X()*nu_inv + xi_[2];    // u = px.X/nu + u0
//            _v = xi_[1]*_P.get_Y()*nu_inv + xi_[3];    // v = py.Y/nu + v0

        }
            break;
          case 1:
            //EUCM
           {
            const double rho = sqrt(xi_[5]*(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y())+_P.get_Z()*_P.get_Z());
            double x =(1-xi_[4])*_P.get_X()/(_P.get_Z()+rho*(xi_[4]/(1-xi_[4])));
            double y =(1-xi_[4])*_P.get_Y()/(_P.get_Z()+rho*(xi_[4]/(1-xi_[4])));

            const double d = x*x+y*y;
            const double gamma1=xi_[0]/(1-xi_[4]);
            const double gamma2=xi_[1]/(1-xi_[4]);
            const double L = 1+K[0]*d+K[1]*d*d+K[4]*d*d*d;
            x = x*L+2*K[2]*x*y+K[3]*(d+2*x*x);
            y = y*L+2*K[2]*(d+2*y*y)+K[3]*x*y;

            _u = gamma1*x + xi_[2];    // u = px.X/nu + u0
            _v = gamma2*y + xi_[3];    // v = py.Y/nu + v0

//            const double nu_inv = 1./(xi_[4]*rho + (1-xi_[4])*_P.get_Z());
//            _u = xi_[0]*_P.get_X()*nu_inv + xi_[2];    // u = px.X/nu + u0
//            _v = xi_[1]*_P.get_Y()*nu_inv + xi_[3];    // v = py.Y/nu + v0
        }
            break;
          case 2:
        {
            const double rho1 = sqrt(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y()+_P.get_Z()*_P.get_Z());
            const double rho2 = sqrt(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y()+(xi_[5]*rho1*_P.get_Z())*(xi_[5]*rho1*_P.get_Z()));
            const double nu_inv=1./(xi_[4]*rho2+(1-xi_[4])*(xi_[5]*rho1+_P.get_Z()));
            double x = (1-xi_[4])*_P.get_X()*nu_inv;
            double y = (1-xi_[4])*_P.get_Y()*nu_inv;

            const double d = x*x+y*y;
            const double gamma1=xi_[0]/(1-xi_[4]);
            const double gamma2=xi_[1]/(1-xi_[4]);
            const double L = 1+K[0]*d+K[1]*d*d+K[4]*d*d*d;
            x = x*L+2*K[2]*x*y+K[3]*(d+2*x*x);
            y = y*L+2*K[2]*(d+2*y*y)+K[3]*x*y;

            _u = gamma1*x + xi_[2];    // u = px.X/nu + u0
            _v = gamma2*y + xi_[3];    // v = py.Y/nu + v0

        }
            break;
          default:
        {std::cout<<"Not exist model"<<std::endl;}

            break;
        }

    }

    // Jacobian corresponding to the intrinsic parameters

    void computeJacobianIntrinsic(const vpPoint &_P, vpMatrix &_J)
    {
        switch (camera_model) {
        case 0:
        {
            const double rho = sqrt(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y()+_P.get_Z()*_P.get_Z());
//            const double rho_inv = 1./rho;
            const double nu_inv = 1./(xi_[4]*rho + (1-xi_[4])*_P.get_Z());
            _J.resize(2,5);
            _J[0][0] = _P.get_X()*nu_inv;                                                                       // du/dpx
            _J[1][1] = _P.get_Y()*nu_inv;                                                                       // dv/dpy
            _J[0][2] = _J[1][3] = 1;                                                                            // du/du0, dv/dv0
            _J[0][4] = xi_[0]*_P.get_X()*(_P.get_Z()-rho)*nu_inv*nu_inv;                                        // du/dalpha
            _J[1][4] = xi_[1]*_P.get_Y()*(_P.get_Z()-rho)*nu_inv*nu_inv;                                        // dv/dalpha
//            _J.cppPrint(std::cout, "Ji");


        }
            break;
        case 1:
        {
            const double rho = sqrt(xi_[5]*(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y())+_P.get_Z()*_P.get_Z());
            const double rho_inv = 1./rho;
            const double nu_inv = 1./(xi_[4]*rho + (1-xi_[4])*_P.get_Z());
            const double nu2rho_inv = rho_inv*nu_inv*nu_inv;

            _J.resize(2,6);
            _J[0][0] = _P.get_X()*nu_inv;                                                                       // du/dpx
            _J[1][1] = _P.get_Y()*nu_inv;                                                                       // dv/dpy
            _J[0][2] = _J[1][3] = 1;                                                                            // du/du0, dv/dv0
            _J[0][4] = xi_[0]*_P.get_X()*(_P.get_Z()-rho)*nu_inv*nu_inv;                                        // du/dalpha
            _J[1][4] = xi_[1]*_P.get_Y()*(_P.get_Z()-rho)*nu_inv*nu_inv;                                        // dv/dalpha
            _J[0][5] = -0.5*xi_[0]*_P.get_X()*xi_[4]*(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y())*nu2rho_inv;  // du/dbeta
            _J[1][5] = -0.5*xi_[1]*_P.get_Y()*xi_[4]*(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y())*nu2rho_inv;  // dv/dbeta
       }
            break;
        case 2:
        {


        }
            break;
        default:
            break;
        }

    }

    void computeJacobianDistortion(const vpPoint &_P, vpMatrix &_J)
    {
        vpMatrix dPdD,dDdV;
        double rho;
        double x,y;
        switch (camera_model) {
        case 0:
           {
            rho = sqrt(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y()+_P.get_Z()*_P.get_Z());
            x = (1-xi_[4])*(1-xi_[4])*_P.get_X()/(_P.get_Z()*(1-xi_[4])+rho*xi_[4]);
            y = (1-xi_[4])*(1-xi_[4])*_P.get_Y()/(_P.get_Z()*(1-xi_[4])+rho*xi_[4]);
        }
            break;
        case 1:
        {
            rho = sqrt(xi_[5]*(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y())+_P.get_Z()*_P.get_Z());
            x = (1-xi_[4])*_P.get_X()/(_P.get_Z()+rho*(xi_[4]/(1-xi_[4])));
            y = (1-xi_[4])*_P.get_Y()/(_P.get_Z()+rho*(xi_[4]/(1-xi_[4])));
        }
            break;
        default:
            break;
        }

        const double d = x*x+y*y;
        dPdD.resize(2,2);
        dDdV.resize(2,5);

        const double gamma1=xi_[0]/(1-xi_[4]);
        const double gamma2=xi_[1]/(1-xi_[4]);
        dPdD[0][0] = gamma1;
        dPdD[1][1] = gamma2;

        dDdV[0][0]=x*d;
        dDdV[0][1]=x*d*d;
        dDdV[0][2]=2*x*y;
        dDdV[0][3]=d+2*x*x;
        dDdV[0][4]=x*d*d*d;

        dDdV[1][0]=y*d;
        dDdV[1][1]=y*d*d;
        dDdV[1][2]=d+2*y*y;
        dDdV[1][3]=2*x*y;
        dDdV[1][4]=y*d*d*d;
        //            dPdD.cppPrint(std::cout, "dPdD");
        //            dDdV.cppPrint(std::cout, "dDdV");
        //            _J.cppPrint(std::cout, "Di");
        _J=dPdD*dDdV;

    }


    // Jacobian wrt extrinsic parameters
    // J should be 2x6
    void computeJacobianExtrinsic(const vpPoint &_P, vpMatrix &_J)
    {
        vpMatrix dPdD,dDdH,dHdX_;
        double rho;
        double x,y;
        const double gamma1=xi_[0]/(1-xi_[4]);
        const double gamma2=xi_[1]/(1-xi_[4]);
        dPdD.resize(2,2);
        dDdH.resize(2,2);
        dHdX_.resize(2,3);

        dPdD[0][0] = gamma1;
        dPdD[1][1] = gamma2;

        switch (camera_model) {
        case 0:
          {
           rho = sqrt(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y()+_P.get_Z()*_P.get_Z());

           x = (1-xi_[4])*(1-xi_[4])*_P.get_X()/(_P.get_Z()*(1-xi_[4])+rho*xi_[4]);
           y = (1-xi_[4])*(1-xi_[4])*_P.get_Y()/(_P.get_Z()*(1-xi_[4])+rho*xi_[4]);



           const double xi=xi_[4]/(1-xi_[4]);
           const double nu_inv = 1./(rho*(_P.get_Z()+xi*rho)*(_P.get_Z()+xi*rho));
           dHdX_[0][0] = nu_inv*(rho*_P.get_Z()+xi*(_P.get_Y()*_P.get_Y()+_P.get_Z()*_P.get_Z()));
           dHdX_[0][1] = -nu_inv*xi*_P.get_X()*_P.get_Y();
           dHdX_[0][2] = nu_inv*_P.get_X()*(-xi*_P.get_Z()-rho);
           dHdX_[1][0] = -nu_inv*xi*_P.get_X()*_P.get_Y();
           dHdX_[1][1] = nu_inv*(rho*_P.get_Z()+xi*(_P.get_X()*_P.get_X()+_P.get_Z()*_P.get_Z()));
           dHdX_[1][2] = nu_inv*_P.get_Y()*(-xi*_P.get_Z()-rho);


        }
            break;
        case 1:
        {
            rho = sqrt(xi_[5]*(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y())+_P.get_Z()*_P.get_Z());

            x =(1-xi_[4])* _P.get_X()/(_P.get_Z()+rho*(xi_[4]/(1-xi_[4])));
            y =(1-xi_[4])* _P.get_Y()/(_P.get_Z()+rho*(xi_[4]/(1-xi_[4])));



            const double xi=xi_[4]/(1-xi_[4]);
            const double nu_inv = 1./(rho*(_P.get_Z()+xi*rho)*(_P.get_Z()+xi*rho));
            dHdX_[0][0] = nu_inv*(rho*_P.get_Z()+xi*(rho*rho-xi_[5]*_P.get_X()*_P.get_X()));
            dHdX_[0][1] = -nu_inv*xi*_P.get_X()*_P.get_Y()*xi_[5];
            dHdX_[0][2] = nu_inv*_P.get_X()*(-xi*_P.get_Z()-rho);
            dHdX_[1][0] = -nu_inv*xi*_P.get_X()*_P.get_Y()*xi_[5];
            dHdX_[1][1] = nu_inv*(rho*_P.get_Z()+xi*(rho*rho-xi_[5]*_P.get_Y()*_P.get_Y()));
            dHdX_[1][2] = nu_inv*_P.get_Y()*(-xi*_P.get_Z()-rho);

        }
            break;
        default:
            break;
        }

        const double d = x*x+y*y;
        const double I=2*K[0]*x*y+4*K[1]*d*x*y+2*K[2]*x+2*K[3]*y+6*K[4]*x*y*d*d;
        dDdH[0][0] = 1+K[0]*(d+2*x*x)+K[1]*d*(d+4*x*x)+2*K[2]*y+6*K[3]*x+K[5]*d*d*(d+6*x*x);
        dDdH[0][1] = dDdH[1][0] = I;
        dDdH[1][1] = 1+K[0]*(d+2*y*y)+K[1]*d*(d+4*y*y)+6*K[2]*y+2*K[3]*x+K[5]*d*d*(d+6*y*y);


        dPdX_=dPdD*dDdH*dHdX_;
//      dPdD.cppPrint(std::cout, "dPdD");
//        dDdH.cppPrint(std::cout, "dDdH");
//      dHdX_.cppPrint(std::cout, "dHdX_");
        vpFeatureBuilder::create(p_, _P);

        _J = dPdX_ * p_.interaction();
//           dPdX_.cppPrint(std::cout, "dPdX_");
//           _J.cppPrint(std::cout, "Li");


//        const double rho = sqrt(xi_[5]*(_P.get_X()*_P.get_X()+_P.get_Y()*_P.get_Y())+_P.get_Z()*_P.get_Z());
//        const double rho_inv = 1./rho;
//        const double nu_inv = 1./(xi_[4]*rho + (1-xi_[4])*_P.get_Z());
//        const double nu2rho_inv = rho_inv*nu_inv*nu_inv;

//        dPdX_[0][0] = xi_[0]*(nu_inv - xi_[4]*xi_[5]*_P.get_X()*_P.get_X()*nu2rho_inv);     // du/dX
//        dPdX_[0][1] = -xi_[0]*(xi_[4]*xi_[5]*_P.get_X()*_P.get_Y()*nu2rho_inv);             // du/dY
//        dPdX_[0][2] = -xi_[0]*_P.get_X()*nu_inv*nu_inv*(1-xi_[4]+xi_[4]*_P.get_Z()*rho_inv);// du/dZ
//        dPdX_[1][0] = -xi_[1]*(xi_[4]*xi_[5]*_P.get_X()*_P.get_Y()*nu2rho_inv);             // dv/dX
//        dPdX_[1][1] = xi_[1]*(nu_inv - xi_[4]*xi_[5]*_P.get_Y()*_P.get_Y()*nu2rho_inv);     // dv/dY
//        dPdX_[1][2] = -xi_[1]*_P.get_Y()*nu_inv*nu_inv*(1-xi_[4]+xi_[4]*_P.get_Z()*rho_inv);// dv/dZ

//        vpFeatureBuilder::create(p_, _P);

//        _J = dPdX_ * p_.interaction();


    }

    // update parameter value
    void updateIntrinsic(const vpColVector &_dxi, const vpColVector &_dk)
    {
        // update
        xi_ += _dxi;
//        K.cppPrint(std::cout, "K1");
        K += _dk;
        _dk.cppPrint(std::cout, "K2");


        // all parameters should be positive
        for(int i = 0;i<6;++i)
        {
            if(xi_[i] < 0)
               { xi_[i] = 0;}
        }
        for(int i = 0;i<5;++i)
        {
            if(K[i] < 0)
                K[i] = 0;
        }
//        K.cppPrint(std::cout, "K");

        // alpha should be lesser than 1
        if(xi_[4] > 1)
        {xi_[4] = 0.99;}

//        xi_.cppPrint(std::cout, "xi");
    }
};

}

#endif // DISTORTION_CAMERA_H


