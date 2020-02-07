#include <visp/vpHomogeneousMatrix.h>

#include <visp/vpPoint.h>
#include <visp/vpSubColVector.h>
#include <visp/vpSubMatrix.h>
#include <visp/vpFeaturePoint.h>
#include <visp/vpFeatureBuilder.h>
#include <visp/vpExponentialMap.h>
#include <visp/vpAdaptiveGain.h>
#include <visp/vpIoTools.h>
#include <fstream>

#include <opencv2/calib3d/calib3d.hpp>

#include <vvs.h>
#include <grid_tracker.h>
#include <perspective_camera.h>
#include <distortion_camera.h>
#include <cb_tracker.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::stringstream;
using cv::waitKey;
using namespace covis;

int main()
{
    // load calibration images from hard drive
    const string base = "/home/zhx/Documents/calibration/images/";

//    const string base = "/home/zhx/Documents/RicohTheta_ws/src/ricoh_camera/calibration/camodocal/back_images/";
    const string prefix = "back000";
//    const string prefix = "fishgrid";
    // init empty vector of detected patterns
    vector<Pattern> patterns;
    patterns.clear();
    patterns.reserve(36);
    bool findCorners=0;

//    GridTracker tracker;      // this tracker detects a 6x6 grid of points
    CBTracker tracker(8,5); // this one is to be given the chessboard dimension (8x6)

    // read images while the corresponding file exists
    // images are displayed to ensure the detection was performed
    int image_num=0;
    while(true)
    {
        stringstream ss;
//        ss << prefix << patterns.size() << ".png";
        ss << prefix << image_num << ".png";
        image_num+=1;
        std::ifstream testfile(base + ss.str());
        if(testfile.good())
        {
            testfile.close();
            Pattern pat;
            pat.im =  cv::imread(base + ss.str());

            findCorners=tracker.detect(pat.im, pat.point);
            pat.window = ss.str();
//            cout<<q<<endl;
            if(findCorners)
            {
                drawSeq(pat.window, pat.im, pat.point);
                patterns.push_back(pat);

            }
            waitKey(1);
            // draw extraction results

        }
        else
            break;
    }
    cout << "Found " << patterns.size() << " images" << endl;

    // create a camera model (Perspective or Distortion)
    // default parameters should be guessed from image dimensions

    const double pxy = 0.5*( patterns.data()->im.rows +patterns.data()->im.cols );
    const double alpha = 0.5;
    const double beta = 0.5;
    DistortionCamera cam(pxy , pxy , 0.5* patterns.data()->im.cols , 0.5*patterns.data()->im.rows, alpha, beta);   // not a very good guess

    // initiate virtual visual servoing with inter-point distance and pattern dimensions
    VVS vvs(cam, 0.03, 8, 5);

//    // calibrate from all images
    vvs.calibrate(patterns);

//    // print results
    cout << "Final calibration: " << cam.xi_.t() << endl;

    // this will wait for a key pressed to stop the program
    waitKey(0);
}
