#pragma once
#include <opencv2/core/core.hpp>    
#include <opencv2/highgui/highgui.hpp> 
#include  <cstring> 
#include <opencv2/imgproc/imgproc.hpp>  
using namespace cv;

class FindLocation
{
public:
	FindLocation();
	~FindLocation();
	int mainTextDetection(Mat byteQueryimage,int co);
	int co;
public:
	
};
namespace TextDection
{
	//struct SWTPoint2d {
	//	int x;
	//	int y;
	//	float SWT;
	//};
	//struct Ray {
	//	SWTPoint2d p;
	//	SWTPoint2d q;
	//	std::vector<SWTPoint2d> points;
	//};
	//bool Point2dSort(SWTPoint2d const & lhs, SWTPoint2d const & rhs) { return lhs.SWT < rhs.SWT; }
	//typedef std::pair<SWTPoint2d, SWTPoint2d> SWTPointPair2d;
	//typedef std::pair<cv::Point, cv::Point>   SWTPointPair2i;
	//struct Point2dFloat {
	//	float x;
	//	float y;
	//};

	//struct Chain {
	//	int p;
	//	int q;
	//	float dist;
	//	bool merged;
	//	Point2dFloat direction;
	//	std::vector<int> components;
	//};

	//struct Point3dFloat {
	//	float x;
	//	float y;
	//	float z;
	//};

	//const Scalar BLUE(255, 0, 0);
	//const Scalar GREEN(0, 255, 0);
	//const Scalar RED(0, 0, 255);
	//bool chainSortDist(const Chain &lhs, const Chain &rhs) 
	//{
	//	return lhs.dist < rhs.dist;
	//}
	//bool chainSortLength(const Chain &lhs, const Chain &rhs)
	//{
	//	return lhs.components.size() > rhs.components.size();
	//}




	//std::vector< std::vector<SWTPoint2d> > findLegallyConnectedComponents(cv::Mat& SWTImage, std::vector<Ray> & rays);

	//std::vector< std::vector<SWTPoint2d> >
	//	findLegallyConnectedComponentsRAY(IplImage * SWTImage,
	//		std::vector<Ray> & rays);

	//void componentStats(IplImage * SWTImage,
	//	const std::vector<SWTPoint2d> & component,
	//	float & mean, float & variance, float & median,
	//	int & minx, int & miny, int & maxx, int & maxy);

	//void filterComponents(cv::Mat& SWTImage,
	//	std::vector<std::vector<SWTPoint2d> > & components,
	//	std::vector<std::vector<SWTPoint2d> > & validComponents,
	//	std::vector<Point2dFloat> & compCenters,
	//	std::vector<float> & compMedians,
	//	std::vector<SWTPoint2d> & compDimensions,
	//	std::vector<SWTPointPair2d > & compBB);

	//std::vector<Chain> makeChains(const cv::Mat& colorImage,
	//	std::vector<std::vector<SWTPoint2d> > & components,
	//	std::vector<Point2dFloat> & compCenters,
	//	std::vector<float> & compMedians,
	//	std::vector<SWTPoint2d> & compDimensions,
	//	std::vector<SWTPointPair2d > & compBB);
}

