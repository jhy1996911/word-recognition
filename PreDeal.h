#pragma once
#include <opencv2/core/core.hpp>    
#include <opencv2/highgui/highgui.hpp> 
#include  <cstring> 
#include <opencv2/imgproc/imgproc.hpp>  
using namespace cv;
class PreDeal
{
public:
	Mat im;
	PreDeal();
	~PreDeal();
	Mat getMat();
	Mat setMat(Mat im);
	void adaptiveThreshold(Mat input, Mat & output, int width, int height);
};
