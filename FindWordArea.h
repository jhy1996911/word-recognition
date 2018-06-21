#pragma once
#include  "opencv2/text.hpp"
#include  "opencv2/highgui.hpp"
#include  "opencv2/imgproc.hpp"
#include  <vector>
#include  <iostream>
#include  <iomanip>
using namespace cv;
using namespace cv::text;

class FindWordArea
{
public:
	FindWordArea();
	~FindWordArea();
	void FindArea(Mat src, std::vector<Rect>&a);

};

