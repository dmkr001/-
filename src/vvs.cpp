
#include <vvs.h>
#include <visp/vpSubColVector.h>
#include <visp/vpExponentialMap.h>

using namespace covis;
using std::cout;
using std::endl;
using std::vector;
using std::string;

/* Calibrate the camera from sequence of images and extracted points
    _pat is a list of image patterns with:
        - _pat[i].im is the i-eth image
        - _pat[i].point are the pixel points extracted from this image
        - _pat[i].window is the name of the window to display the results

*/
void VVS::calibrate(std::vector<Pattern> &_pat)
{
    // number of images
    const int n = _pat.size();
    // number of points
    const int m = r_*c_;
    // number of parameters
    const int n_k=5;
    const int n_xi = cam_->nbParam();
    // the decrease speed
    const double lamda=0.3;
    // the initial guess for M
    vector<vpHomogeneousMatrix> M(n);
    for(unsigned int k=0;k<n;++k)
        M[k].buildFrom(
            0,0.,0.5,                                                                                   // translation
            0,0,atan2(_pat[k].point[5].y-_pat[k].point[0].y, _pat[k].point[5].x-_pat[k].point[0].x));   // rotation

    // current and desired features: 2 coord. for each m points for each n images
    vpColVector s(2*m*n), sd(2*m*n);

    // write measured positions as desired features
    unsigned int row = 0;

    unsigned int k = 0;                     // image number
    for(auto &pat: _pat)                    // loop through all images
    {
        for(unsigned int i=0;i<m;++i)       // loop through all points
        {
            row = 2*i+2*m*k;                // corresponding row in the global error vector
            sd[row] = pat.point[i].x;       // x value
            sd[row+1] =  pat.point[i].y;    // y value
        }
        k++;
    }

    // update vector

    vpColVector dx(n_xi + n_k + 6*n);      // intrinsic parameters + 6 velocities for each image
//    vpColVector dx(n_xi + 6*n);
//    vpSubColVector dxi(dx, 0, n_xi); // part of dx that represents the intrinsic parameters
//    vpSubColVector dk(dx, n_xi, n_k);
    // Jacobians
    vpMatrix J(2*m*n, n_xi+n_k+6*n);     // global Jacobian
//    vpMatrix J(2*m*n, n_xi+6*n);
    vpMatrix Ji(2, n_xi);            // Jacobian for 1 point for intrinsic
    vpMatrix Li(2, 6);               // Jacobian for 1 point for extrinsic
    vpMatrix Di(2, n_k);
    unsigned int iter=0;
    const int camera_model=0;

    dx = 1; // set to 1 so that it is higher than minimum error at first
    while(dx.euclideanNorm() > 0.00001 && iter++ < 101)
    {
        /* first we have to compute for all points from all images:
         * - the pixel coordinates corresponding to the current estimation of intrinsic and extrinsic
         * - the intrinsic Jacobian
         * - the extrinsic Jacobian
         */

        unsigned int k = 0;                     // image number
        for(auto &pat: _pat)                    // loop through all images with index k
        {
            for(unsigned int i=0;i<m;++i)       // loop through all points with index i
            {
                // corresponding row in the global error vector
                row = 2*i+2*m*k;
                // track the point X_[i] from the current estimation of M_k
                X_[i].track(M[k]);
                // use the current intrinsic estimation to project into pixels
                cam_->project(X_[i],s[row],s[row+1]);

                // compute the intrinsic Jacobian for this point
                cam_->computeJacobianIntrinsic(X_[i],Ji);

                cam_->computeJacobianDistortion(X_[i],Di);

                // compute the extrinsic Jacobian for this point
                cam_->computeJacobianExtrinsic(X_[i],Li);

                // write Ji and Li inside J, using the putAt function
                putAt(J, Ji, 2*i+2*m*k, 0);   // writes Ji in J
                putAt(J, Di, 2*i+2*m*k, n_xi);
                putAt(J, Li, 2*i+2*m*k, n_xi+n_k+k*6);   // write Li in J

//                putAt(J, Li, 2*i+2*m*k, n_xi+k*6);
            }
            k++;
        }


        /* Here we have the error (s-sd) and the corresponding Jacobian J
         * We can compute an update for the unknown with the pseudo inverse
         */

        dx = -lamda*J.pseudoInverse()*(s-sd);
//        dx.cppPrint(std::cout, "dx");


        // we now update the intrinsic parameters, common to all images
        vpSubColVector dxi(dx, 0, n_xi); // part of dx that represents the intrinsic parameters
        vpSubColVector dk(dx, n_xi, n_k);
        cam_->updateIntrinsic(dxi, dk);
//        dxi.cppPrint(std::cout, "dxi");
//        dk.cppPrint(std::cout, "dk:");



        // and the extrinsic parameters for each image
        for(unsigned int k=0;k<n;++k)
        {
            vpSubColVector v(dx, n_xi+n_k+6*k, 6);          // extract velocity vector from global update dx
//            vpSubColVector v(dx, n_xi+6*k, 6);
//            v.cppPrint(std::cout, "v");
            M[k] = vpExponentialMap::direct(v).inverse() * M[k];
        }

        // if we want to see what's happening during the loop
        display(_pat, M, s, 5);

    }

    // here the camera should be calibrated and all poses M_k should be estimated, display the result
    display(_pat, M, s, 0);
    vpColVector err=(s-sd);
    std::cout<<err.frobeniusNorm()/err.getRows()<<endl;

}


/* Compute the pose of the camera from a single image and extracted points and a guess on the current pose M
    _pat is an image patterns with:
        - _pat.im is the image
        - _pat.point are the pixel points extracted from this image
        - _pat.window is the name of the window to display the results

    _reset can be used as an external way to tell the VVS not to use the last M found but to reinitialize one from a wild guess
            (used when the pose was badly calibrated)
*/

void VVS::computePose(Pattern &_pat, vpHomogeneousMatrix &_M, const bool &_reset)
{







}
