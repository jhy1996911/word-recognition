#include "stdafx.h"
#include "PreDeal.h"
#include "OCRProject.h"
#include "OCRProjectDlg.h"
#include "afxdialogex.h"
using namespace cv;
PreDeal::PreDeal()
{
}
PreDeal::~PreDeal()
{
}


Mat PreDeal::getMat()
{
	return im;
}

Mat PreDeal::setMat(Mat img)
{
	im = img;
}

void PreDeal::adaptiveThreshold(Mat input, Mat &output, int width, int height)
{
	if (input.channels() != 1)
	{
		cv::cvtColor(input, input, CV_BGR2GRAY);
	}
	int S = width >> 3;
	int T = 15;
	unsigned long* integralImg = 0;
	int i, j;
	long sum = 0;
	int count = 0;
	int index;
	int x1, y1, x2, y2;
	int s2 = S / 2;

	output.create(input.size(), input.type());
	// create the integral image  
	integralImg = (unsigned long*)malloc(width*height * sizeof(unsigned long*));

	for (i = 0; i<width; i++)
	{
		// reset this column sum  
		sum = 0;
		const uchar* inData = input.ptr<uchar>(i);
		for (j = 0; j<height; j++)
		{
			index = j*width + i;
			sum += inData[j];
			if (i == 0)
				integralImg[index] = sum;
			else
				integralImg[index] = integralImg[index - 1] + sum;
		}
	}
	for (i = 0; i<width; i++)
	{
		const uchar* inData = input.ptr<uchar>(i);
		uchar* indata = output.ptr<uchar>(i);
		for (j = 0; j<height; j++)
		{
			index = j*width + i;
			// set the SxS region  
			x1 = i - s2; x2 = i + s2;
			y1 = j - s2; y2 = j + s2;
			// check the border  
			if (x1 < 0) x1 = 0;
			if (x2 >= width) x2 = width - 1;
			if (y1 < 0) y1 = 0;
			if (y2 >= height) y2 = height - 1;
			count = (x2 - x1)*(y2 - y1);
			sum = integralImg[y2*width + x2] -
				integralImg[y1*width + x2] -
				integralImg[y2*width + x1] +
				integralImg[y1*width + x1];
			if ((long)((inData[j])* count) < (long)(sum*(100 - T) / 100))
				indata[j] = 0;
			else
				indata[j] = 255;
		}
	}
	imshow("1", output);
}
