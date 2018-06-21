#include "stdafx.h"
#include "FindLocation.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/unordered_map.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <math.h>
#include <time.h>
#include <utility>
#include <algorithm>
#include <vector>
using namespace TextDection;
#define PI 3.14159265
FindLocation::FindLocation()
{
}


FindLocation::~FindLocation()
{
}

struct SWTPoint2d {
	int x;
	int y;
	float SWT;
};
struct Ray {
	SWTPoint2d p;
	SWTPoint2d q;
	std::vector<SWTPoint2d> points;
};
bool Point2dSort(SWTPoint2d const & lhs, SWTPoint2d const & rhs) { return lhs.SWT < rhs.SWT; }
typedef std::pair<SWTPoint2d, SWTPoint2d> SWTPointPair2d;
typedef std::pair<cv::Point, cv::Point>   SWTPointPair2i;
struct Point2dFloat {
	float x;
	float y;
};

struct Chain {
	int p;
	int q;
	float dist;
	bool merged;
	Point2dFloat direction;
	std::vector<int> components;
};

struct Point3dFloat {
	float x;
	float y;
	float z;
};

const Scalar BLUE(255, 0, 0);
const Scalar GREEN(0, 255, 0);
const Scalar RED(0, 0, 255);
bool chainSortDist(const Chain &lhs, const Chain &rhs)
{
	return lhs.dist < rhs.dist;
}
bool chainSortLength(const Chain &lhs, const Chain &rhs)
{
	return lhs.components.size() > rhs.components.size();
}
void strokeWidthTransform(const Mat& edgeImage,Mat& gradientX,Mat& gradientY,bool dark_on_light,Mat& SWTImage,std::vector<Ray> & rays) {
	// First pass
	float prec = .05;
	//遍历像素点，如果该点的为边缘，则将该点的坐标保存下来
	for (int row = 0; row < edgeImage.rows; row++) {
		const uchar* ptr = (const uchar*)edgeImage.ptr(row);
		for (int col = 0; col < edgeImage.cols; col++) {
			if (*ptr > 0) {   
				Ray r;
				SWTPoint2d p;
				p.x = col;
				p.y = row;
				r.p = p;
				std::vector<SWTPoint2d> points;
				points.push_back(p);

				float curX = (float)col + 0.5;
				float curY = (float)row + 0.5;
				int curPixX = col;
				int curPixY = row;
				float G_x = gradientX.at<float>(row, col);
				float G_y = gradientY.at<float>(row, col);
				// 标准化梯度
				float mag = sqrt((G_x * G_x) + (G_y * G_y));
				if (dark_on_light) {
					G_x = -G_x / mag;
					G_y = -G_y / mag;
				}
				else {
					G_x = G_x / mag;
					G_y = G_y / mag;

				}
				while (true) {
					curX += G_x*prec;
					curY += G_y*prec;
					if ((int)(floor(curX)) != curPixX || (int)(floor(curY)) != curPixY) {
						curPixX = (int)(floor(curX));
						curPixY = (int)(floor(curY));
						//判断每个点是否在图像外
						if (curPixX < 0 || (curPixX >= SWTImage.cols) || curPixY < 0 || (curPixY >= SWTImage.rows)) {
							break;
						}
						SWTPoint2d pnew;
						pnew.x = curPixX;
						pnew.y = curPixY;
						points.push_back(pnew);

						//二次删选
						if (edgeImage.at<uchar>(curPixY, curPixX) > 0) {
							r.q = pnew;
							// dot product
							float G_xt = gradientX.at<float>(curPixY, curPixX);
							float G_yt = gradientY.at<float>(curPixY, curPixX);
							//梯度标准化
							mag = sqrt((G_xt * G_xt) + (G_yt * G_yt));
							if (dark_on_light) {
								G_xt = -G_xt / mag;
								G_yt = -G_yt / mag;
							}
							else {
								G_xt = G_xt / mag;
								G_yt = G_yt / mag;

							}
							//寻找梯度方向反方向的点，如果梯度方向加上反方向的绝对值小于90度，则认定为笔画
							if (acos(G_x * -G_xt + G_y * -G_yt) < PI / 2.0) {
								//计算笔画宽度
								float length = sqrt(((float)r.q.x - (float)r.p.x)*((float)r.q.x - (float)r.p.x) + ((float)r.q.y - (float)r.p.y)*((float)r.q.y - (float)r.p.y));
								//遍历SWTPoint2d中所有的点，将路径上所有的点都赋值为该距离
								for (std::vector<SWTPoint2d>::iterator pit = points.begin(); pit != points.end(); pit++) 
								{
									if (SWTImage.at<float>(pit->y, pit->x) < 0) {
										SWTImage.at<float>(pit->y, pit->x) = length;
									}
									else {
										//如果路径上的点的距离已经存在，则取最小距离
										SWTImage.at<float>(pit->y, pit->x) = std::min(length, SWTImage.at<float>(pit->y, pit->x));
									}
								}
								//该点存入r
								r.points = points;
								rays.push_back(r);
							}
							break;
						}
					}
				}
			}
			ptr++;
		}
	}

}


//第二次扫描
//对于所有第一次扫描到的路径，
//求出路径上的中值，将所有这条路径上大于中值的点全部赋值为中值。

void SWTMedianFilter(Mat& SWTImage, std::vector<Ray> & rays) {
	for (auto& rit : rays) {
		for (auto& pit : rit.points) {
			pit.SWT = SWTImage.at<float>(pit.y, pit.x);
		}
		std::sort(rit.points.begin(), rit.points.end(), &Point2dSort);
		float median = (rit.points[rit.points.size() / 2]).SWT;
		for (auto& pit : rit.points) {
			SWTImage.at<float>(pit.y, pit.x) = std::min(pit.SWT, median);
		}
	}
}

//图像标准化
void normalizeImage(const Mat& input, Mat& output) {
	assert(input.depth() == CV_32F);
	assert(input.channels() == 1);
	assert(output.depth() == CV_32F);
	assert(output.channels() == 1);

	float maxVal = 0;
	float minVal = 1e100;
	for (int row = 0; row < input.rows; row++) {
		const float* ptr = (const float*)input.ptr(row);
		for (int col = 0; col < input.cols; col++) {
			if (*ptr < 0) {}
			else {
				maxVal = std::max(*ptr, maxVal);
				minVal = std::min(*ptr, minVal);
			}
			ptr++;
		}
	}

	float difference = maxVal - minVal;
	for (int row = 0; row < input.rows; row++) {
		const float* ptrin = (const float*)input.ptr(row);
		float* ptrout = (float*)output.ptr(row);
		for (int col = 0; col < input.cols; col++) {
			if (*ptrin < 0) {
				*ptrout = 1;
			}
			else {
				*ptrout = ((*ptrin) - minVal) / difference;
			}
			ptrout++;
			ptrin++;
		}
	}
}


std::vector< std::vector<SWTPoint2d> >
findLegallyConnectedComponentsRAY(Mat& SWTImage,
	std::vector<Ray> & rays)
{
	boost::unordered_map<int, int> map;
	boost::unordered_map<int, SWTPoint2d> revmap;

	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
	int num_vertices = 0;
	// Number vertices for graph.  Associate each point with number
	for (int row = 0; row < SWTImage.rows; row++) {
		float * ptr = (float*)SWTImage.ptr(row);
		for (int col = 0; col < SWTImage.cols; col++) {
			if (*ptr > 0) {
				map[row * SWTImage.cols + col] = num_vertices;
				SWTPoint2d p;
				p.x = col;
				p.y = row;
				revmap[num_vertices] = p;
				num_vertices++;
			}
			ptr++;
		}
	}

	Graph g(num_vertices);

	// Traverse and add edges to graph
	for (std::vector<Ray>::const_iterator it = rays.begin(); it != rays.end(); it++) {
		float lastSW = 0;
		int lastRow = 0;
		int lastCol = 0;
		for (std::vector<SWTPoint2d>::const_iterator it2 = it->points.begin(); it2 != it->points.end(); it2++) {
			float currentSW = SWTImage.at<float>(it2->y, it2->x);
			if (lastSW == 0) {}
			else if (lastSW / currentSW <= 3.0 || currentSW / lastSW <= 3.0) {
				boost::add_edge(map.at(it2->y * SWTImage.cols + it2->x), map.at(lastRow * SWTImage.cols + lastCol), g);
			}
			lastSW = currentSW;
			lastRow = it2->y;
			lastCol = it2->x;
		}
		lastSW = 0;
		lastRow = 0;
		lastCol = 0;
	}

	std::vector<int> c(num_vertices);

	int num_comp = connected_components(g, &c[0]);

	std::vector<std::vector<SWTPoint2d> > components;
	components.reserve(num_comp);
	std::cout << "Before filtering, " << num_comp << " components and " << num_vertices << " vertices" << std::endl;
	for (int j = 0; j < num_comp; j++) {
		std::vector<SWTPoint2d> tmp;
		components.push_back(tmp);
	}
	for (int j = 0; j < num_vertices; j++) {
		SWTPoint2d p = revmap[j];
		(components[c[j]]).push_back(p);
	}

	return components;
}

//计算笔画宽度的均值，方差，中值
void componentStats(Mat& SWTImage,
	const std::vector<SWTPoint2d> & component,
	float & mean, float & variance, float & median,
	int & minx, int & miny, int & maxx, int & maxy)
{
	std::vector<float> temp;
	temp.reserve(component.size());
	mean = 0;
	variance = 0;
	minx = 1000000;
	miny = 1000000;
	maxx = 0;
	maxy = 0;
	for (std::vector<SWTPoint2d>::const_iterator it = component.begin(); it != component.end(); it++) {
		float t = SWTImage.at<float>(it->y, it->x);
		mean += t;
		temp.push_back(t);
		miny = std::min(miny, it->y);
		minx = std::min(minx, it->x);
		maxy = std::max(maxy, it->y);
		maxx = std::max(maxx, it->x);
	}
	mean = mean / ((float)component.size());
	for (std::vector<float>::const_iterator it = temp.begin(); it != temp.end(); it++) {
		variance += (*it - mean) * (*it - mean);
	}
	variance = variance / ((float)component.size());
	std::sort(temp.begin(), temp.end());
	median = temp[temp.size() / 2];
}


void filterComponents(Mat& SWTImage,
	std::vector<std::vector<SWTPoint2d> > & components,
	std::vector<std::vector<SWTPoint2d> > & validComponents,
	std::vector<Point2dFloat> & compCenters,
	std::vector<float> & compMedians,
	std::vector<SWTPoint2d> & compDimensions,
	std::vector<SWTPointPair2d > & compBB)
{
	validComponents.reserve(components.size());
	compCenters.reserve(components.size());
	compMedians.reserve(components.size());
	compDimensions.reserve(components.size());
	// bounding boxes
	compBB.reserve(components.size());
	for (std::vector<std::vector<SWTPoint2d> >::iterator it = components.begin(); it != components.end(); it++) {
		// 计算笔画宽度均值, 方差, 中值
		float mean, variance, median;
		int minx, miny, maxx, maxy;
		//统计三个要素
		componentStats(SWTImage, (*it), mean, variance, median, minx, miny, maxx, maxy);

		// 检查方差是否小于均值
		if (variance > 0.5 * mean) {
			continue;
		}

		float length = (float)(maxx - minx + 1);
		float width = (float)(maxy - miny + 1);

		// 检查字体高度
		if (width > 300) {
			continue;
		}

		float area = length * width;
		float rminx = (float)minx;
		float rmaxx = (float)maxx;
		float rminy = (float)miny;
		float rmaxy = (float)maxy;
		// 计算旋转的边界框
		float increment = 1. / 36.;
		for (float theta = increment * PI; theta<PI / 2.0; theta += increment * PI) {
			float xmin, xmax, ymin, ymax, xtemp, ytemp, ltemp, wtemp;
			xmin = 1000000;
			ymin = 1000000;
			xmax = 0;
			ymax = 0;
			for (unsigned int i = 0; i < (*it).size(); i++) {
				xtemp = (*it)[i].x * cos(theta) + (*it)[i].y * -sin(theta);
				ytemp = (*it)[i].x * sin(theta) + (*it)[i].y * cos(theta);
				xmin = std::min(xtemp, xmin);
				xmax = std::max(xtemp, xmax);
				ymin = std::min(ytemp, ymin);
				ymax = std::max(ytemp, ymax);
			}
			ltemp = xmax - xmin + 1;
			wtemp = ymax - ymin + 1;
			if (ltemp*wtemp < area) {
				area = ltemp*wtemp;
				length = ltemp;
				width = wtemp;
			}
		}
		// 检验连通区域的高宽比是否在 1/10 和 10 之间
		if (length / width < 1. / 10. || length / width > 10.) {
			continue;
		}

		// 计算连通区域的直径
		// 计算连通区域的密度
		std::vector <std::vector<float> > denseRepr;
		denseRepr.reserve(maxx - minx + 1);
		for (int i = 0; i < maxx - minx + 1; i++) {
			std::vector<float> tmp;
			tmp.reserve(maxy - miny + 1);
			denseRepr.push_back(tmp);
			for (int j = 0; j < maxy - miny + 1; j++) {
					denseRepr[i].push_back(0);
			}
		}
		for (std::vector<SWTPoint2d>::iterator pit = it->begin(); pit != it->end(); pit++) {
			(denseRepr[pit->x - minx])[pit->y - miny] = 1;
		}
		// 标记连通区域
		const int num_nodes = it->size();
		/*
		E edges[] = { E(0,2),
		E(1,1), E(1,3), E(1,4),
		E(2,1), E(2,3),
		E(3,4),
		E(4,0), E(4,1) };

		Graph G(edges + sizeof(edges) / sizeof(E), weights, num_nodes);
		*/
		Point2dFloat center;
		center.x = ((float)(maxx + minx)) / 2.0;
		center.y = ((float)(maxy + miny)) / 2.0;

		SWTPoint2d dimensions;
		dimensions.x = maxx - minx + 1;
		dimensions.y = maxy - miny + 1;

		SWTPoint2d bb1;
		bb1.x = minx;
		bb1.y = miny;

		SWTPoint2d bb2;
		bb2.x = maxx;
		bb2.y = maxy;
		SWTPointPair2d pair(bb1, bb2);

		compBB.push_back(pair);
		compDimensions.push_back(dimensions);
		compMedians.push_back(median);
		compCenters.push_back(center);
		validComponents.push_back(*it);
	}
	std::vector<std::vector<SWTPoint2d > > tempComp;
	std::vector<SWTPoint2d > tempDim;
	std::vector<float > tempMed;
	std::vector<Point2dFloat > tempCenters;
	std::vector<SWTPointPair2d > tempBB;
	tempComp.reserve(validComponents.size());
	tempCenters.reserve(validComponents.size());
	tempDim.reserve(validComponents.size());
	tempMed.reserve(validComponents.size());
	tempBB.reserve(validComponents.size());
	for (unsigned int i = 0; i < validComponents.size(); i++) {
		int count = 0;
		for (unsigned int j = 0; j < validComponents.size(); j++) {
			if (i != j) {
				if (compBB[i].first.x <= compCenters[j].x && compBB[i].second.x >= compCenters[j].x &&
					compBB[i].first.y <= compCenters[j].y && compBB[i].second.y >= compCenters[j].y) {
					count++;
				}
			}
		}
		if (count < 2) {
			tempComp.push_back(validComponents[i]);
			tempCenters.push_back(compCenters[i]);
			tempMed.push_back(compMedians[i]);
			tempDim.push_back(compDimensions[i]);
			tempBB.push_back(compBB[i]);
		}
	}
	validComponents = tempComp;
	compDimensions = tempDim;
	compMedians = tempMed;
	compCenters = tempCenters;
	compBB = tempBB;

	compDimensions.reserve(tempComp.size());
	compMedians.reserve(tempComp.size());
	compCenters.reserve(tempComp.size());
	validComponents.reserve(tempComp.size());
	compBB.reserve(tempComp.size());

	std::cout << "After filtering " << validComponents.size() << " components" << std::endl;
}

void renderComponents(const Mat& SWTImage, std::vector<std::vector<SWTPoint2d> > & components, Mat& output) {
	output.setTo(0);

	for (auto& component : components) {
		for (auto& pit : component) {
			output.at<float>(pit.y, pit.x) = SWTImage.at<float>(pit.y, pit.x);
		}
	}
	for (int row = 0; row < output.rows; row++) {
		float* ptr = (float*)output.ptr(row);
		for (int col = 0; col < output.cols; col++) {
			if (*ptr == 0) {
				*ptr = -1;
			}
			ptr++;
		}
	}
	float maxVal = 0;
	float minVal = 1e100;
	for (int row = 0; row < output.rows; row++) {
		const float* ptr = (const float*)output.ptr(row);
		for (int col = 0; col < output.cols; col++) {
			if (*ptr == 0) {}
			else {
				maxVal = std::max(*ptr, maxVal);
				minVal = std::min(*ptr, minVal);
			}
			ptr++;
		}
	}
	float difference = maxVal - minVal;
	for (int row = 0; row < output.rows; row++) {
		float* ptr = (float*)output.ptr(row);
		for (int col = 0; col < output.cols; col++) {
			if (*ptr < 1) {
				*ptr = 1;
			}
			else {
				*ptr = ((*ptr) - minVal) / difference;
			}
			ptr++;
		}
	}

}

//展示联通区域
void renderComponentsWithBoxes(Mat& SWTImage, std::vector<std::vector<SWTPoint2d> > & components,
	std::vector<SWTPointPair2d > & compBB, Mat& output) {
	Mat outTemp(output.size(), CV_32FC1);

	renderComponents(SWTImage, components, outTemp);

	std::vector<SWTPointPair2i> bb;
	bb.reserve(compBB.size());
	for (auto& it : compBB) {
		Point2i p0 = cvPoint(it.first.x, it.first.y);
		Point2i p1 = cvPoint(it.second.x, it.second.y);
		SWTPointPair2i pair(p0, p1);
		bb.push_back(pair);
	}

	Mat out(output.size(), CV_8UC1);
	outTemp.convertTo(out, CV_8UC1, 255.);
	cvtColor(out, output, CV_GRAY2RGB);

	int count = 0;
	for (auto it : bb) {
		Scalar c;
		if (count % 3 == 0) {
			c = BLUE;
		}
		else if (count % 3 == 1) {
			c = GREEN;
		}
		else {
			c = RED;
		}
		count++;
		rectangle(output, it.first, it.second, c, 2);
	}
}

//链接各个连通区域
bool sharesOneEnd(Chain c0, Chain c1) {
	if (c0.p == c1.p || c0.p == c1.q || c0.q == c1.q || c0.q == c1.p) {
		return true;
	}
	else {
		return false;
	}
}
std::vector<Chain> makeChains(const Mat& colorImage,
	std::vector<std::vector<SWTPoint2d> > & components,
	std::vector<Point2dFloat> & compCenters,
	std::vector<float> & compMedians,
	std::vector<SWTPoint2d> & compDimensions,
	std::vector<SWTPointPair2d > & compBB) {
	assert(compCenters.size() == components.size());
	// make vector of color averages
	std::vector<Point3dFloat> colorAverages;
	colorAverages.reserve(components.size());
	for (std::vector<std::vector<SWTPoint2d> >::iterator it = components.begin(); it != components.end(); it++) {
		Point3dFloat mean;
		mean.x = 0;
		mean.y = 0;
		mean.z = 0;
		int num_points = 0;
		for (std::vector<SWTPoint2d>::iterator pit = it->begin(); pit != it->end(); pit++) {
			mean.x += (float)colorImage.at<uchar>(pit->y, (pit->x) * 3);
			mean.y += (float)colorImage.at<uchar>(pit->y, (pit->x) * 3 + 1);
			mean.z += (float)colorImage.at<uchar>(pit->y, (pit->x) * 3 + 2);
			num_points++;
		}
		mean.x = mean.x / ((float)num_points);
		mean.y = mean.y / ((float)num_points);
		mean.z = mean.z / ((float)num_points);
		colorAverages.push_back(mean);
	}

	// form all eligible pairs and calculate the direction of each
	std::vector<Chain> chains;
	for (unsigned int i = 0; i < components.size(); i++) {
		for (unsigned int j = i + 1; j < components.size(); j++) {
			// TODO add color metric
			if ((compMedians[i] / compMedians[j] <= 2.0 || compMedians[j] / compMedians[i] <= 2.0) &&
				(compDimensions[i].y / compDimensions[j].y <= 2.0 || compDimensions[j].y / compDimensions[i].y <= 2.0)) {
				float dist = (compCenters[i].x - compCenters[j].x) * (compCenters[i].x - compCenters[j].x) +
					(compCenters[i].y - compCenters[j].y) * (compCenters[i].y - compCenters[j].y);
				float colorDist = (colorAverages[i].x - colorAverages[j].x) * (colorAverages[i].x - colorAverages[j].x) +
					(colorAverages[i].y - colorAverages[j].y) * (colorAverages[i].y - colorAverages[j].y) +
					(colorAverages[i].z - colorAverages[j].z) * (colorAverages[i].z - colorAverages[j].z);
				if (dist < 9 * (float)(std::max(std::min(compDimensions[i].x, compDimensions[i].y), std::min(compDimensions[j].x, compDimensions[j].y)))
					*(float)(std::max(std::min(compDimensions[i].x, compDimensions[i].y), std::min(compDimensions[j].x, compDimensions[j].y)))
					&& colorDist < 1600) {
					Chain c;
					c.p = i;
					c.q = j;
					std::vector<int> comps;
					comps.push_back(c.p);
					comps.push_back(c.q);
					c.components = comps;
					c.dist = dist;
					float d_x = (compCenters[i].x - compCenters[j].x);
					float d_y = (compCenters[i].y - compCenters[j].y);
					/*
					float d_x = (compBB[i].first.x - compBB[j].second.x);
					float d_y = (compBB[i].second.y - compBB[j].second.y);
					*/
					float mag = sqrt(d_x*d_x + d_y*d_y);
					d_x = d_x / mag;
					d_y = d_y / mag;
					Point2dFloat dir;
					dir.x = d_x;
					dir.y = d_y;
					c.direction = dir;
					chains.push_back(c);

					/*std::cerr << c.p << " " << c.q << std::endl;
					std::cerr << c.direction.x << " " << c.direction.y << std::endl;
					std::cerr << compCenters[c.p].x << " " << compCenters[c.p].y << std::endl;
					std::cerr << compCenters[c.q].x << " " << compCenters[c.q].y << std::endl;
					std::cerr << std::endl;
					std::cerr << colorDist << std::endl; */
				}
			}
		}
	}
	std::cout << chains.size() << " eligible pairs" << std::endl;
	std::sort(chains.begin(), chains.end(), &chainSortDist);

	std::cerr << std::endl;
	const float strictness = PI / 6.0;
	//merge chains
	int merges = 1;
	while (merges > 0) {
		for (unsigned int i = 0; i < chains.size(); i++) {
			chains[i].merged = false;
		}
		merges = 0;
		std::vector<Chain> newchains;
		for (unsigned int i = 0; i < chains.size(); i++) {
			for (unsigned int j = 0; j < chains.size(); j++) {
				if (i != j) {
					if (!chains[i].merged && !chains[j].merged && sharesOneEnd(chains[i], chains[j])) {
						if (chains[i].p == chains[j].p) {
							if (acos(chains[i].direction.x * -chains[j].direction.x + chains[i].direction.y * -chains[j].direction.y) < strictness) {
								/*      if (chains[i].p == chains[i].q || chains[j].p == chains[j].q) {
								std::cout << "CRAZY ERROR" << std::endl;
								} else if (chains[i].p == chains[j].p && chains[i].q == chains[j].q) {
								std::cout << "CRAZY ERROR" << std::endl;
								} else if (chains[i].p == chains[j].q && chains[i].q == chains[j].p) {
								std::cout << "CRAZY ERROR" << std::endl;
								}
								std::cerr << 1 <<std::endl;

								std::cerr << chains[i].p << " " << chains[i].q << std::endl;
								std::cerr << chains[j].p << " " << chains[j].q << std::endl;
								std::cerr << compCenters[chains[i].q].x << " " << compCenters[chains[i].q].y << std::endl;
								std::cerr << compCenters[chains[i].p].x << " " << compCenters[chains[i].p].y << std::endl;
								std::cerr << compCenters[chains[j].q].x << " " << compCenters[chains[j].q].y << std::endl;
								std::cerr << std::endl; */

								chains[i].p = chains[j].q;
								for (std::vector<int>::iterator it = chains[j].components.begin(); it != chains[j].components.end(); it++) {
									chains[i].components.push_back(*it);
								}
								float d_x = (compCenters[chains[i].p].x - compCenters[chains[i].q].x);
								float d_y = (compCenters[chains[i].p].y - compCenters[chains[i].q].y);
								chains[i].dist = d_x * d_x + d_y * d_y;

								float mag = sqrt(d_x*d_x + d_y*d_y);
								d_x = d_x / mag;
								d_y = d_y / mag;
								Point2dFloat dir;
								dir.x = d_x;
								dir.y = d_y;
								chains[i].direction = dir;
								chains[j].merged = true;
								merges++;
								/*j=-1;
								i=0;
								if (i == chains.size() - 1) i=-1;
								std::stable_sort(chains.begin(), chains.end(), &chainSortLength);*/
							}
						}
						else if (chains[i].p == chains[j].q) {
							if (acos(chains[i].direction.x * chains[j].direction.x + chains[i].direction.y * chains[j].direction.y) < strictness) {
								/*
								if (chains[i].p == chains[i].q || chains[j].p == chains[j].q) {
								std::cout << "CRAZY ERROR" << std::endl;
								} else if (chains[i].p == chains[j].p && chains[i].q == chains[j].q) {
								std::cout << "CRAZY ERROR" << std::endl;
								} else if (chains[i].p == chains[j].q && chains[i].q == chains[j].p) {
								std::cout << "CRAZY ERROR" << std::endl;
								}
								std::cerr << 2 <<std::endl;

								std::cerr << chains[i].p << " " << chains[i].q << std::endl;
								std::cerr << chains[j].p << " " << chains[j].q << std::endl;
								std::cerr << chains[i].direction.x << " " << chains[i].direction.y << std::endl;
								std::cerr << chains[j].direction.x << " " << chains[j].direction.y << std::endl;
								std::cerr << compCenters[chains[i].q].x << " " << compCenters[chains[i].q].y << std::endl;
								std::cerr << compCenters[chains[i].p].x << " " << compCenters[chains[i].p].y << std::endl;
								std::cerr << compCenters[chains[j].p].x << " " << compCenters[chains[j].p].y << std::endl;
								std::cerr << std::endl; */

								chains[i].p = chains[j].p;
								for (std::vector<int>::iterator it = chains[j].components.begin(); it != chains[j].components.end(); it++) {
									chains[i].components.push_back(*it);
								}
								float d_x = (compCenters[chains[i].p].x - compCenters[chains[i].q].x);
								float d_y = (compCenters[chains[i].p].y - compCenters[chains[i].q].y);
								float mag = sqrt(d_x*d_x + d_y*d_y);
								chains[i].dist = d_x * d_x + d_y * d_y;

								d_x = d_x / mag;
								d_y = d_y / mag;

								Point2dFloat dir;
								dir.x = d_x;
								dir.y = d_y;
								chains[i].direction = dir;
								chains[j].merged = true;
								merges++;
								/*j=-1;
								i=0;
								if (i == chains.size() - 1) i=-1;
								std::stable_sort(chains.begin(), chains.end(), &chainSortLength); */
							}
						}
						else if (chains[i].q == chains[j].p) {
							if (acos(chains[i].direction.x * chains[j].direction.x + chains[i].direction.y * chains[j].direction.y) < strictness) {
								/*                           if (chains[i].p == chains[i].q || chains[j].p == chains[j].q) {
								std::cout << "CRAZY ERROR" << std::endl;
								} else if (chains[i].p == chains[j].p && chains[i].q == chains[j].q) {
								std::cout << "CRAZY ERROR" << std::endl;
								} else if (chains[i].p == chains[j].q && chains[i].q == chains[j].p) {
								std::cout << "CRAZY ERROR" << std::endl;
								}
								std::cerr << 3 <<std::endl;

								std::cerr << chains[i].p << " " << chains[i].q << std::endl;
								std::cerr << chains[j].p << " " << chains[j].q << std::endl;

								std::cerr << compCenters[chains[i].p].x << " " << compCenters[chains[i].p].y << std::endl;
								std::cerr << compCenters[chains[i].q].x << " " << compCenters[chains[i].q].y << std::endl;
								std::cerr << compCenters[chains[j].q].x << " " << compCenters[chains[j].q].y << std::endl;
								std::cerr << std::endl; */
								chains[i].q = chains[j].q;
								for (std::vector<int>::iterator it = chains[j].components.begin(); it != chains[j].components.end(); it++) {
									chains[i].components.push_back(*it);
								}
								float d_x = (compCenters[chains[i].p].x - compCenters[chains[i].q].x);
								float d_y = (compCenters[chains[i].p].y - compCenters[chains[i].q].y);
								float mag = sqrt(d_x*d_x + d_y*d_y);
								chains[i].dist = d_x * d_x + d_y * d_y;


								d_x = d_x / mag;
								d_y = d_y / mag;
								Point2dFloat dir;
								dir.x = d_x;
								dir.y = d_y;

								chains[i].direction = dir;
								chains[j].merged = true;
								merges++;
								/*j=-1;
								i=0;
								if (i == chains.size() - 1) i=-1;
								std::stable_sort(chains.begin(), chains.end(), &chainSortLength); */
							}
						}
						else if (chains[i].q == chains[j].q) {
							if (acos(chains[i].direction.x * -chains[j].direction.x + chains[i].direction.y * -chains[j].direction.y) < strictness) {
								/*           if (chains[i].p == chains[i].q || chains[j].p == chains[j].q) {
								std::cout << "CRAZY ERROR" << std::endl;
								} else if (chains[i].p == chains[j].p && chains[i].q == chains[j].q) {
								std::cout << "CRAZY ERROR" << std::endl;
								} else if (chains[i].p == chains[j].q && chains[i].q == chains[j].p) {
								std::cout << "CRAZY ERROR" << std::endl;
								}
								std::cerr << 4 <<std::endl;
								std::cerr << chains[i].p << " " << chains[i].q << std::endl;
								std::cerr << chains[j].p << " " << chains[j].q << std::endl;
								std::cerr << compCenters[chains[i].p].x << " " << compCenters[chains[i].p].y << std::endl;
								std::cerr << compCenters[chains[i].q].x << " " << compCenters[chains[i].q].y << std::endl;
								std::cerr << compCenters[chains[j].p].x << " " << compCenters[chains[j].p].y << std::endl;
								std::cerr << std::endl; */
								chains[i].q = chains[j].p;
								for (std::vector<int>::iterator it = chains[j].components.begin(); it != chains[j].components.end(); it++) {
									chains[i].components.push_back(*it);
								}
								float d_x = (compCenters[chains[i].p].x - compCenters[chains[i].q].x);
								float d_y = (compCenters[chains[i].p].y - compCenters[chains[i].q].y);
								chains[i].dist = d_x * d_x + d_y * d_y;

								float mag = sqrt(d_x*d_x + d_y*d_y);
								d_x = d_x / mag;
								d_y = d_y / mag;
								Point2dFloat dir;
								dir.x = d_x;
								dir.y = d_y;
								chains[i].direction = dir;
								chains[j].merged = true;
								merges++;
								/*j=-1;
								i=0;
								if (i == chains.size() - 1) i=-1;
								std::stable_sort(chains.begin(), chains.end(), &chainSortLength);*/
							}
						}
					}
				}
			}
		}
		for (unsigned int i = 0; i < chains.size(); i++) {
			if (!chains[i].merged) {
				newchains.push_back(chains[i]);
			}
		}
		chains = newchains;
		std::stable_sort(chains.begin(), chains.end(), &chainSortLength);
	}

	std::vector<Chain> newchains;
	newchains.reserve(chains.size());
	for (std::vector<Chain>::iterator cit = chains.begin(); cit != chains.end(); cit++) {
		if (cit->components.size() >= 3) {
			newchains.push_back(*cit);
		}
	}
	chains = newchains;
	std::cout << chains.size() << " chains after merging" << std::endl;
	return chains;
}

//展示字符块连通区域
void renderChains(Mat& SWTImage,
	std::vector<std::vector<SWTPoint2d> > & components,
	std::vector<Chain> & chains,
	Mat& output) {
	// keep track of included components
	std::vector<bool> included;
	included.reserve(components.size());
	for (unsigned int i = 0; i != components.size(); i++) {
		included.push_back(false);
	}
	for (std::vector<Chain>::iterator it = chains.begin(); it != chains.end(); it++) {
		for (std::vector<int>::iterator cit = it->components.begin(); cit != it->components.end(); cit++) {
			included[*cit] = true;
		}
	}
	std::vector<std::vector<SWTPoint2d> > componentsRed;
	for (unsigned int i = 0; i != components.size(); i++) {
		if (included[i]) {
			componentsRed.push_back(components[i]);
		}
	}
	std::cout << componentsRed.size() << " components after chaining" << std::endl;
	Mat outTemp(output.size(), CV_32FC1);
	renderComponents(SWTImage, componentsRed, outTemp);
	outTemp.convertTo(output, CV_8UC1, 255);

}
std::vector< std::vector<SWTPoint2d> > findLegallyConnectedComponents(Mat& SWTImage, std::vector<Ray> & rays) {
	boost::unordered_map<int, int> map;
	boost::unordered_map<int, SWTPoint2d> revmap;

	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
	int num_vertices = 0;
	// Number vertices for graph.  Associate each point with number
	for (int row = 0; row < SWTImage.rows; row++) {
		float * ptr = (float*)SWTImage.ptr(row);
		for (int col = 0; col < SWTImage.cols; col++) {
			if (*ptr > 0) {
				map[row * SWTImage.cols + col] = num_vertices;
				SWTPoint2d p;
				p.x = col;
				p.y = row;
				revmap[num_vertices] = p;
				num_vertices++;
			}
			ptr++;
		}
	}

	Graph g(num_vertices);

	for (int row = 0; row < SWTImage.rows; row++) {
		float * ptr = (float*)SWTImage.ptr(row);
		for (int col = 0; col < SWTImage.cols; col++) {
			if (*ptr > 0) {
				// check pixel to the right, right-down, down, left-down
				int this_pixel = map[row * SWTImage.cols + col];
				if (col + 1 < SWTImage.cols) {
					float right = SWTImage.at<float>(row, col + 1);
					if (right > 0 && ((*ptr) / right <= 3.0 || right / (*ptr) <= 3.0))
						boost::add_edge(this_pixel, map.at(row * SWTImage.cols + col + 1), g);
				}
				if (row + 1 < SWTImage.rows) {
					if (col + 1 < SWTImage.cols) {
						float right_down = SWTImage.at<float>(row + 1, col + 1);
						if (right_down > 0 && ((*ptr) / right_down <= 3.0 || right_down / (*ptr) <= 3.0))
							boost::add_edge(this_pixel, map.at((row + 1) * SWTImage.cols + col + 1), g);
					}
					float down = SWTImage.at<float>(row + 1, col);
					if (down > 0 && ((*ptr) / down <= 3.0 || down / (*ptr) <= 3.0))
						boost::add_edge(this_pixel, map.at((row + 1) * SWTImage.cols + col), g);
					if (col - 1 >= 0) {
						float left_down = SWTImage.at<float>(row + 1, col - 1);
						if (left_down > 0 && ((*ptr) / left_down <= 3.0 || left_down / (*ptr) <= 3.0))
							boost::add_edge(this_pixel, map.at((row + 1) * SWTImage.cols + col - 1), g);
					}
				}
			}
			ptr++;
		}
	}

	std::vector<int> c(num_vertices);

	int num_comp = connected_components(g, &c[0]);

	std::vector<std::vector<SWTPoint2d> > components;
	components.reserve(num_comp);
	std::cout << "Before filtering, " << num_comp << " components and " << num_vertices << " vertices" << std::endl;
	for (int j = 0; j < num_comp; j++) {
		std::vector<SWTPoint2d> tmp;
		components.push_back(tmp);
	}
	for (int j = 0; j < num_vertices; j++) {
		SWTPoint2d p = revmap[j];
		(components[c[j]]).push_back(p);
	}

	return components;
}




Mat textDetection(const Mat& input, bool dark_on_light) {
	assert(input.depth() == CV_8U);
	assert(input.channels() == 3);
	// 灰度化处理
	Mat grayimage(input.size(), CV_8UC1);
//	cvtColor(input, grayimage, CV_RGB2GRAY);
	grayimage = input.clone();
	// Create Canny image
	double threshold_low = 175;
	double threshold_high = 320;
	Mat edgeimage(input.size(), CV_8UC1);
	Canny(grayimage, edgeimage, threshold_low, threshold_high, 3); //边缘检测
	imwrite("canny.jpg", edgeimage);
	// Create gradient X, gradient Y
	Mat gaussianimage(input.size(), CV_32FC1);
	grayimage.convertTo(gaussianimage, CV_32FC1, 1. / 255.);
	GaussianBlur(gaussianimage, gaussianimage, Size(5, 5), 0);
	Mat gradientX(input.size(), CV_32FC1);
	Mat gradientY(input.size(), CV_32FC1);
	Scharr(gaussianimage, gradientX, -1, 1, 0);
	Scharr(gaussianimage, gradientY, -1, 0, 1);
	GaussianBlur(gradientX, gradientX, Size(3, 3), 0);
	GaussianBlur(gradientY, gradientY, Size(3, 3), 0);

	// Calculate SWT and return ray vectors
	std::vector<Ray> rays;
	Mat SWTimage(input.size(), CV_32FC1);
	for (int row = 0; row < input.rows; row++) {
		float* ptr = (float*)SWTimage.ptr(row);
		for (int col = 0; col < input.cols; col++) {
			*ptr++ = -1;
		}
	}//初始化SWT图像

	strokeWidthTransform(edgeimage, gradientX, gradientY, dark_on_light, SWTimage, rays);
	SWTMedianFilter(SWTimage, rays);

	Mat output2(input.size(), CV_32FC1);
	normalizeImage(SWTimage, output2);
	Mat saveSWT(input.size(), CV_8UC1);
	output2.convertTo(saveSWT, CV_8UC1, 255);
	imwrite("SWT.png", saveSWT);



	//计算合理的联通要素在swt图像和梯度图像
	//返回容器的容器，外部容器表示连通区域，内部容器表示联通区域的点
	std::vector<std::vector<SWTPoint2d> > components = findLegallyConnectedComponents(SWTimage, rays);
	// 过滤掉不相关区域
	std::vector<std::vector<SWTPoint2d> > validComponents;
	std::vector<SWTPointPair2d > compBB;
	std::vector<Point2dFloat> compCenters;
	std::vector<float> compMedians;
	std::vector<SWTPoint2d> compDimensions;
	filterComponents(SWTimage, components, validComponents, compCenters, compMedians, compDimensions, compBB);

	Mat output3(input.size(), CV_8UC3);
	//展示联通区域
	renderComponentsWithBoxes(SWTimage, validComponents, compBB, output3);
	imwrite("components.png", output3);
	//

	// 链接连通区域
	std::vector<Chain> chains;
	chains = makeChains(input, validComponents, compCenters, compMedians, compDimensions, compBB);

	Mat output4(input.size(), CV_8UC1);
	//展示字符块连通区域
	renderChains(SWTimage, validComponents, chains, output4);
	imwrite ( "text.png", output4);
	//imshow("1", output4);
	//waitKey(0);
	Mat output5(input.size(), CV_8UC3);
	cvtColor(output4, output5, CV_GRAY2RGB);


	/*Iplimage * output =
	cvCreateimage ( input.size(), CV_8UC3 );
	renderChainsWithBoxes ( SWTimage, validComponents, chains, compBB, output); */
	return output5;
}


int FindLocation:: mainTextDetection(Mat byteQueryimage,int co) {
	if (byteQueryimage.empty() ||  co == -1) {
		return -1;
	}
	// Detect text in the image
	imshow("1", byteQueryimage);
	imwrite("src.jpg", byteQueryimage);
	Mat output = textDetection(byteQueryimage, co);
	imwrite("result.jpg", output);
	imshow("2", output);
	waitKey(0);
	return 0;
}