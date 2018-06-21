#include "stdafx.h"
#include "FindWordArea.h"
using namespace std;

FindWordArea::FindWordArea()
{
}


FindWordArea::~FindWordArea()
{
}


bool Judgment_intersection(Rect a, Rect b)
{
	Rect c = a&b;
	if (c.width == 0 && c.height == 0 && c.x == 0 && c.y == 0)
	{
		return false;
	}
	else
		return true;
}


void groups_draw(Mat &src, vector<Rect> &groups)
{
	vector<Rect> groups_res;
	int flag;
	int flag1;
	while (1)
	{
		flag = 0;
		for (int i = 0; i < groups.size(); ++i)
		{
			flag1 = 0;
			for (int j = i + 1; j < groups.size(); ++j)
			{
				if (Judgment_intersection(groups[i], groups[j]))
				{
					Rect c;
					//c.x = min((double)groups[i].x,(double)groups[j].x);
					//c.y = min((double)groups[i].y, (double)groups[j].y);
					//c.width = max(groups[i].x + groups[i].width, groups[j].x + groups[j].width)-c.x;
					//c.height=max(groups[i].y + groups[i].height, groups[j].y + groups[j].height)-c.y;
					c = groups[i] | groups[j];
					if (i < j)
					{
						groups.erase(groups.begin() + i);
						groups.erase(groups.begin() + j - 1);
					}
					else
					{
						groups.erase(groups.begin() + j);
						groups.erase(groups.begin() + i - 1);
					}
					groups.push_back(c);
					flag = 1;
					flag1 = 1;
					break;
				}
				if (flag1)
					break;
			}
		}
		if (!flag)
		{
			break;
		}
	}
	for (int i = 0; i < groups.size(); ++i)
	{
		if (groups[i].x + groups[i].width>src. cols)
		{
			groups[i].x = src.cols- groups[i].width;
		}
		else if (groups[i].y + groups[i].height>src.rows)
		{
			groups[i].y = src.rows - groups[i].height;
		}
	}
	groups[0].x = 310;
	groups[0].width = 280;
	for (int i = (int)groups.size() - 1; i >= 0; i--)
	{
		if (src.type() == CV_8UC3)
			rectangle(src, groups.at(i).tl(), groups.at(i).br(), Scalar(0, 255, 255), 3, 8);
		else
			rectangle(src, groups.at(i).tl(), groups.at(i).br(), Scalar(255), 3, 8);
	}
}



void er_show(vector<Mat> &channels, vector<vector<ERStat> > &regions)
{
	for (int c = 0; c<(int)channels.size(); c++)
	{
		Mat dst = Mat::zeros(channels[0].rows + 2, channels[0].cols + 2, CV_8UC1);
		for (int r = 0; r<(int)regions[c].size(); r++)
		{
			ERStat er = regions[c][r];
			if (er.parent != NULL) // deprecate the root region
			{
				int newMaskVal = 255;
				int flags = 4 + (newMaskVal << 8) + FLOODFILL_FIXED_RANGE + FLOODFILL_MASK_ONLY;
				floodFill(channels[c], dst, Point(er.pixel%channels[c].cols, er.pixel / channels[c].cols),
					Scalar(255), 0, Scalar(er.level), Scalar(0), flags);
			}
		}
		char buff[10]; char *buff_ptr = buff;
		sprintf_s(buff, "channel %d", c);
		imshow(buff_ptr, dst);
	}
	waitKey(-1);
}
void FindWordArea::FindArea(Mat src, vector<Rect>& r)
{
	// Extract channels to be processed individually
	vector<Mat> channels;
	computeNMChannels(src, channels);

	int cn = (int)channels.size();
	// Append negative channels to detect ER- (bright regions over dark background)
	for (int c = 0; c < cn - 1; c++)
		channels.push_back(255 - channels[c]);

	// Create ERFilter objects with the 1st and 2nd stage default classifiers
	Ptr<ERFilter> er_filter1 = createERFilterNM1(loadClassifierNM1("trained_classifierNM1.xml"), 16, 0.00015f, 0.13f, 0.2f, true, 0.1f);
	Ptr<ERFilter> er_filter2 = createERFilterNM2(loadClassifierNM2("trained_classifierNM2.xml"), 0.5);

	vector<vector<ERStat> > regions(channels.size());
	// Apply the default cascade classifier to each independent channel (could be done in parallel)
	cout << "Extracting Class Specific Extremal Regions from " << (int)channels.size() << " channels ..." << endl;
	cout << "    (...) this may take a while (...)" << endl << endl;
	for (int c = 0; c<(int)channels.size(); c++)
	{
		er_filter1->run(channels[c], regions[c]);
		er_filter2->run(channels[c], regions[c]);
	}

	// Detect character groups
	cout << "Grouping extracted ERs ... ";
	vector< vector<Vec2i> > region_groups;
	vector<Rect> groups_boxes;
	erGrouping(src, channels, regions, region_groups, groups_boxes, ERGROUPING_ORIENTATION_HORIZ);
	//erGrouping(src, channels, regions, region_groups, groups_boxes, ERGROUPING_ORIENTATION_ANY, "./trained_classifier_erGrouping.xml", 0.5);

	// draw groups
	groups_draw(src, groups_boxes);
	r.swap(groups_boxes);
	imshow("grouping", src);

	cout << "Done!" << endl << endl;
	cout << "Press 'space' to show the extracted Extremal Regions, any other key to exit." << endl << endl;
	if ((waitKey() & 0xff) == ' ')
		er_show(channels, regions);

	// memory clean-up
	er_filter1.release();
	er_filter2.release();
	regions.clear();
	if (!groups_boxes.empty())
	{
		groups_boxes.clear();
	}
}
void show_help_and_exit(const char *cmd)
{
	cout << "    Usage: " << cmd << " <input_image> " << endl;
	cout << "    Default classifier files (trained_classifierNM*.xml) must be in current directory" << endl << endl;
	exit(-1);
}
bool panduan_xiangjiao(Rect a, Rect b)
{
	// 需要排除特殊情况：一个矩形在另一个矩形内  
	if ((a.y < b.y && a.y - a.height > b.y - b.height && a.x + a.width < b.x + b.width && a.x > b.x) || (b.y < a.y && b.y - b.height > a.y - a.height && b.x + b.width < a.x + a.width && b.x > a.x))
	{
		return false;
	}
	if ((a.y< b.y - b.height && a.x + a.width < b.x) || (a.y - a.height > b.y && a.x > b.x + b.width))
	{
		return false;
	}
	return true;
}
