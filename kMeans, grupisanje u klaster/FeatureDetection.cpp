#include "FeatureDetection.h"
#include "SIFTLib.h"
#include "ColorSpaces.h"
#include "ImageFilter.h"
#include <list>

using namespace std;

void SIFTDetect(uchar input[], int xSize, int ySize)
{
	SiftKeypointList kptList;

	uchar* Y_buff = new uchar[xSize*ySize];
	char* U_buff = new char[xSize*ySize / 4];
	char* V_buff = new char[xSize*ySize / 4];

	/* Convert input image to YUV420 image */
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	/* Perform SIFT transformation  */
	kptList = calculateSIFT(Y_buff, xSize, ySize);
	
	/* Mark all keypoints in input image */
	markSIFTKeypointsRGB(input, xSize, ySize, kptList, 255, 0, 0);

	delete[] Y_buff;
	delete[] U_buff;
	delete[] V_buff;
}

void splitFeatures(const SiftKeypointList& in, int width, SiftKeypointList& leftImageKP, SiftKeypointList& rightImageKP)
{
	for (SiftKeypoint kp : in)
	{
		if (kp.c < width / 2)
			leftImageKP.push_back(kp);
		else
			rightImageKP.push_back(kp);
	}
}

double l2Distance(SiftKeypoint kp1, SiftKeypoint kp2)
{
	double ret = 0;
	for (int i = 0; i < DEGREE_OF_DESCRIPTORS; i++)
	{
		ret += (kp1.descriptors[i] - kp2.descriptors[i]) * (kp1.descriptors[i] - kp2.descriptors[i]);
	}
	
	ret = sqrt(ret);
	return ret;
}

void matchFeatures(SiftKeypointList leftImageKP, SiftKeypointList rightImageKP, list<pair<QPoint, QPoint>>& matchPairs, double threshold)
{
	for (SiftKeypoint kp1 : leftImageKP)
	{
		for (SiftKeypoint kp2 : rightImageKP)
		{
			if (l2Distance(kp1, kp2) < threshold)
			{
				matchPairs.push_back(pair<QPoint, QPoint>(QPoint(kp1.r, kp1.c), QPoint(kp2.r, kp2.c)));
			}
		}
	}
}

void SIFTDetectPlusMatch(uchar input[], int xSize, int ySize, double threshold)
{
	SiftKeypointList kptList, kptListLeft, kptListRight;

	uchar* Y_buff = new uchar[xSize*ySize];
	char* U_buff = new char[xSize*ySize / 4];
	char* V_buff = new char[xSize*ySize / 4];

	/* Convert input image to YUV420 image */
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	/* Perform SIFT transformation  */
	kptList = calculateSIFT(Y_buff, xSize, ySize);
	
	/* Split the features of left and right images in separate lists */
	splitFeatures(kptList, xSize, kptListLeft, kptListRight);

	/* Match the features from two images */
	list<pair<QPoint, QPoint>> matchedDots;
	matchFeatures(kptListLeft, kptListRight, matchedDots, threshold);

	/* Draw a line for each mached feature pair */
	drawMatchedFeaturesLinesRGB(input, xSize, ySize, matchedDots, 0, 0, 255);

	delete[] Y_buff;
	delete[] U_buff;
	delete[] V_buff;

}

bool isLocalMaximum(double input[], int xSize, int ySize, int x, int y, int N)
{
	//force odd number
	if (N % 2 == 0)
	{
		N += 1;
	}
	for (int i = x - N / 2; i <= x + N / 2; i++)
	{
		for (int j = y - N / 2; j <= y + N / 2; j++)
		{
			if (input[j*xSize + i] > input[y*xSize + x])
				return false;
		}
	}
	return true;
}


void HarisCornerDetection(uchar input[], int xSize, int ySize, double threshold)
{
	SiftKeypointList kptList, kptListLeft, kptListRight;
	const double alpha = 0.05;
	uchar* sobelVertical = new uchar[xSize*ySize];
	uchar* sobelHorizontal = new uchar[xSize*ySize];
	char* U_buff = new char[xSize*ySize / 4];
	char* V_buff = new char[xSize*ySize / 4];

	/* Convert input image to YUV420 image */
	RGBtoYUV420(input, xSize, ySize, sobelVertical, U_buff, V_buff);
	memcpy(sobelHorizontal, sobelVertical, xSize*ySize);

	double sobelv[9] = { -0.25, -0.5, -0.25, 0, 0, 0, 0.25, 0.5, 0.25 };
	double sobelh[9] = { -0.25, 0, 0.25, -0.5, 0, 0.5, -0.25, 0, 0.25 };

	convolve2D(sobelVertical, xSize, ySize, sobelv, 3);
	convolve2D(sobelHorizontal, xSize, ySize, sobelh, 3);
	
	double* R = new double[xSize*ySize];
	for(int i = 0; i < xSize; i++)
	{
		for(int j = 0; j < ySize; j++)
		{
				double M00 = sobelHorizontal[j*xSize+ i] * sobelHorizontal[j*xSize+ i];
				double M01 = sobelVertical[j*xSize+ i] * sobelHorizontal[j*xSize+ i];
				double M10 = sobelVertical[j*xSize+ i] * sobelHorizontal[j*xSize+ i];
				double M11 = sobelVertical[j*xSize+ i] * sobelVertical[j*xSize+ i];
				
				R[j*xSize+ i] = M00*M11 - M01*M10 - alpha*(M00+M11)*(M00+M11);
		}
	}

	
/* .... */

	delete[] R;
	delete[] sobelVertical;
	delete[] sobelHorizontal;
	delete[] U_buff;
	delete[] V_buff;

}


