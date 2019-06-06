#include "FeatureDetection.h"
#include "SIFTLib.h"
#include "ColorSpaces.h"
#include "ImageFilter.h"
#include <list>

using namespace std;

/*******************************************************************************************************************************/
/* Find SIFT keypoints and mark each one with red color*/
//prepoznaje objekte na slici i boji ih razlicitim bojama
/*******************************************************************************************************************************/
void SIFTDetect(uchar input[], int xSize, int ySize)
{
	SiftKeypointList kptList;
	
	/* Convert input image to YUV420 image */
	uchar* Y_buff = new uchar[xSize * ySize];
	char* U_buff = new char[xSize * ySize / 4];
	char* V_buff = new char[xSize * ySize / 4];
	RGBtoYUV420(input, xSize, ySize, Y_buff, U_buff, V_buff);

	/* Perform SIFT transformation  */
	kptList = calculateSIFT(Y_buff, xSize, ySize);

	/* Mark all keypoints in input image */
	uchar R, G = 255, B;
	markSIFTKeypointsRGB(input, xSize, ySize, kptList, R, G, B);

	delete[] Y_buff;
	delete[] U_buff;
	delete[] V_buff;
}

/*******************************************************************************************************************************/
/* Helper function. Splits image feature set in half (used for "two image" input)  
/*******************************************************************************************************************************/
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

/*******************************************************************************************************************************/
/* Calculate Euclidean between two descriptor vectors
/*******************************************************************************************************************************/
double l2Distance(SiftKeypoint kp1, SiftKeypoint kp2)
{
	double dist = 0.0;
	for (int i = 0; i < DEGREE_OF_DESCRIPTORS; i++)
	{
		dist += pow((kp1.descriptors[i] - kp2.descriptors[i]), 2);
	}
	dist = sqrt(dist);
	return dist;
}

/*******************************************************************************************************************************/
/* Go through the first set of keypoints, and for each keypoint try to find the match in the second set
/* One feature matches the other if Euclidean distance between them is lower than threshold. 
/*******************************************************************************************************************************/
void matchFeatures(SiftKeypointList leftImageKP, SiftKeypointList rightImageKP, list<pair<QPoint, QPoint>>& matchPairs, double threshold)
{
	for (auto point1 : leftImageKP)
	{
		for (auto point2 : rightImageKP)
		{
			if (l2Distance(point1, point2) < threshold)
			{
				matchPairs.push_back(make_pair(QPoint(point1.r, point1.r), QPoint(point2.c, point2.r)));
			}
		}
	}
}

/*******************************************************************************************************************************/
/* Find SIFT keypoints, split the image in half, match the features from one image with those from another
/* one and connect them with green lines
//podeli sliku SplitFeature
//pronadje zajednicke tacke, matchFeatures
//i crta ih
/*******************************************************************************************************************************/
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
	uchar R, G=255, B;
	drawMatchedFeaturesLinesRGB(input, xSize, ySize, matchedDots, R, G, B);

	delete[] Y_buff;
	delete[] U_buff;
	delete[] V_buff;

}

// Basic algorithm for image interpolation.
//interpolacije slike, proslediti zeljenu velicinu
static void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	double scaleX = (double)newXSize / xSize;
	double scaleY = (double)newYSize / ySize;
	int roundX, roundY;
	for (int p = 0; p < newXSize; p++)
	{
		for (int q = 0; q < newYSize; q++)
		{
			roundX = round((double)p / scaleX);
			roundY = round((double)q / scaleY);
			for (int k = 0; k < 3; k++)
			{
				// RGB
				output[q * 3 * newXSize + p * 3 + k] = input[roundY * 3 * xSize + roundX * 3 + k];
			}
		}
	}
}

// Basic algorithm for image rotation.
static void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	int xprim, yprim;

	angle = angle * (3.14 / 180);
	double sinAngle = sin(angle);
	double cosAngle = cos(angle);

	double mfactorSin = m * sinAngle;
	double mfactorCos = m * cosAngle;

	double nfactorSin = n * sinAngle;
	double nfactorCos = n * cosAngle;

	for (int x = 0; x < xSize; x++)
	{
		for (int y = 0; y < ySize; y++)
		{
			xprim = round(x*cosAngle - y * sinAngle - mfactorCos + nfactorSin + m);
			yprim = round(y*cosAngle + x * sinAngle - mfactorSin - nfactorCos + n);

			if (xprim < 0 || xprim >= xSize || yprim < 0 || yprim >= ySize)
			{
				// Black pixel.
				for (int k = 0; k < 3; k++)
				{
					output[y * 3 * xSize + x * 3 + k] = 0;
				}
			}
			else
			{
				// Sample and hold.
				for (int k = 0; k < 3; k++)
				{
					output[y * 3 * xSize + x * 3 + k] = input[yprim * 3 * xSize + xprim * 3 + k];
				}
			}
		}
	}
}

/*
	Calculate difference in size and orientation between left and right image,
	and scale and rotate right image so that it matches the left image.
*/
void equalizeImages(uchar input[], int xSize, int ySize, double threshold)
{
	// TODO: Size difference.
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

	// Pick minimum 5 random pairs and calculate average distance between
	// sift keypoints from left and right image.
	qsrand(0);

	int pairsSize = matchedDots.size();
	// Average distance by x and y coordinate.
	// y - row, 0 index
	// x - columnm 1-index
	double distanceLeft[2] = { 0.0, 0.0}, distanceRight[2] = { 0.0 , 0.0};
	// Ratio between left and right average distance x and y component.
	// y-ratio | 0-index,
	// x-ratio | 1-index
	double ratio[2]; 

	// Helper variables.
	int selected[5] = { -1, -1, -1, -1, -1 }; // remember array for index's of already selected pairs.
	int counter = 0;
	int yprev[2]; // y coordinates from previous iteration, 0-left and 1-right.
	int xprev[2]; // x coordinates from previous iteration, 0-left and 1-right.
	while (counter < 5)
	{
		// Random position of pair.
		int index = qrand() % pairsSize;
		// Check if pair was already selected.
		bool found = false;
		for (int i = 0; i < 5; i++)
		{
			if (selected[i] != -1 && selected[i] == index)
			{
				found = true;
				break;
			}
		}
		
		if (!found)
		{
			// Store index in remember array.
			for (int i = 0; i < 5; i++)
			{
				if (selected[i] != -1)
					selected[i] = index;
			}
			// Forward iterator to selected pair.
			auto it = matchedDots.begin();
			for (int i = 0; i < index; i++)
				it++;

			// Calculate distances.
			if (counter > 1)
			{
				distanceLeft[0] = (distanceLeft[0] + abs(yprev[0] - it->first.y)) / 2; // y-left
				distanceLeft[1] = (distanceLeft[1] + abs(xprev[0] - it->first.x)) / 2; // x-left
				distanceRight[0] = (distanceRight[0] + abs(yprev[1] - it->second.y)) / 2; // y-right
				distanceRight[1] = (distanceRight[1] + abs(xprev[1] - it->second.x)) / 2; // x-right
			}

			yprev[0] = it->first.y; 
			yprev[1] = it->second.y;
			xprev[0] = it->first.x; 
			xprev[1] = it->second.x;
			counter++;
		}
	}

	// Calculate ratio.
	ratio[0] = (distanceLeft[0] / distanceRight[0]); // y-ratio.
	ratio[1] = (distanceLeft[1] / distanceRight[1]); // x-ratio.


	// TODO: Scale right image.
	// Interpolate or decimate right image based on ratio.
	int newYSize = ySize * ratio[0];
	int newXSize = (xSize / 2) * ratio[1];
	int xHalf = xSize / 2;

	// Copy right image.
	uchar* rightImageCopy = new uchar[ySize * xHalf];
	for (int i = xHalf; i < xSize; i++)
	{
		for (int j = 0; j < ySize; j++)
		{
			rightImageCopy[j*xHalf + (i - xHalf)] = input[j*xSize + i];
		}
	}

	// Buffer used for scaling.
	uchar* imageBuffer = new uchar[newXSize * newYSize];

	// Perform image scaling.
	sampleAndHold(rightImageCopy, xHalf, ySize, imageBuffer, newXSize, newYSize);

	// Copy image from buffer to right-image copy and extend borders if needed.
	int yDiff = ySize - newYSize;
	int xDiff = xHalf - newXSize;

	if (yDiff > 0 && xDiff > 0)
	{
		// Extended borders because image buffer is smaller that right image.
		// User larger difference to extend.
		int maxDiff = (yDiff > xDiff) ? yDiff : xDiff;
		uchar* extendedBuffer = new uchar[(newXSize + maxDiff) * (newYSize + maxDiff)];
		extendBorders(imageBuffer, newXSize, newYSize, extendedBuffer, ceil(maxDiff / 2));
		for (int i = 0; i < xHalf; i++)
		{
			for (int j = 0; j < ySize; j++)
			{
				rightImageCopy[j*xHalf + i] = extendedBuffer[j*xHalf + j];
			}
		}
		delete[] extendedBuffer;
	}
	else
	{
		for (int i = 0; i < xHalf; i++)
		{
			for (int j = 0; j < ySize; j++)
			{
				rightImageCopy[j*xHalf + i] = imageBuffer[j*xHalf + j];
			}
		}
	}

	// Place copy of right image into original image.
	for (int i = xHalf; i < xSize; i++)
	{
		for (int j = 0; j < ySize; j++)
		{
			input[j*xSize + i] = rightImageCopy[j*xHalf + (i - xHalf)];
		}
	}


	// TODO: Orientation difference.
	// TODO: Rotate right image.

	delete[] imageBuffer;
	delete[] rightImageCopy;
	delete[] Y_buff;
	delete[] U_buff;
	delete[] V_buff;
}

/*******************************************************************************************************************************/
/* Harris algorithm for corner detection 
//pronalazi uglove
/*******************************************************************************************************************************/
void HarisCornerDetection(uchar input[], int xSize, int ySize, double threshold)
{
	uchar* sobelVertical = new uchar[xSize*ySize];
	uchar* sobelHorizontal = new uchar[xSize*ySize];
	char* U_buff = new char[xSize*ySize / 4];
	char* V_buff = new char[xSize*ySize / 4];

	/* Convert input image to YUV420 image */
	RGBtoYUV420(input, xSize, ySize, sobelVertical, U_buff, V_buff);
	
	/* Create a copy of Y component, since it is needed to calculate derivative in both axis */
	memcpy(sobelHorizontal, sobelVertical, xSize*ySize);

	/* Filter both images with corresponding Sobel operator */
	double sobelVerticalFilter[9] = {-0.25, -0.5, -0.25, 0, 0, 0, 0.25, 0.5, 0.25};
	double sobelHorizontalFilter[9] = { -0.25, 0, 0.25, -0.5, 0, 0.5, -0.25, 0, 0.25};
	convolve2D(sobelHorizontal, xSize, ySize, sobelHorizontalFilter, 9);
	convolve2D(sobelVertical, xSize, ySize, sobelVerticalFilter, 9);

	/* For each pixel, calculate the matrix M, then calculate the R factor and place it in new matrix */
	/* Constant alpha is 0.05. */
	double alpha = 0.05;
	double M[4];
	double* R = new double[xSize * ySize];

	for (int i = 0; i < xSize; i++)
	{
		for (int j = 0; j < ySize; j++)
		{
			M[0] = sobelHorizontal[i * ySize + j] * sobelHorizontal[i * ySize + j];
			M[1] = sobelHorizontal[i * ySize + j] * sobelVertical[i * ySize + j];
			M[2] = sobelHorizontal[i * ySize + j] * sobelVertical[i * ySize + j];
			M[3] = sobelVertical[i * ySize + j] * sobelVertical[i * ySize + j];

			R[i*ySize + j] = (M[0] * M[3] - M[1] * M[2]) - alpha * pow((M[0] + M[3]), 2);
		}
	}

	/* For each entry in R matrix, if the value is greater then threshold, check the 3x3 block arround the pixel
	/* and if it is local maximum, colour the entire 3x3 blok in the input image in blue */
	for (int i = 0; i < xSize; i++)
	{
		for (int j = 0; j < ySize; j++)
		{
			if (R[i*ySize + j] > threshold)
			{
				bool foundLarger = false;
				for (int xx = -1; xx < 2 && !foundLarger; xx++)
				{
					for (int yy = -1; yy < 2; yy++)
					{
						if (i + xx < 0 || i + xx >= xSize)
							break; // Out of bounds.
						else if (j + yy >= 0 && j + yy < ySize)
						{
							// Inside matrix.
							if (R[(i+xx)*ySize + (j+yy)] > R[i*ySize + j])
							{
								// Not local maximum.
								foundLarger = true;
								break;
							}
						}
					}
				}

				if (!foundLarger)
				{
					// Local maximum.
					for (int xx = -1; xx < 2; xx++)
					{
						for (int yy = -1; yy < 2; yy++)
						{
							if (i + xx < 0 || i + xx >= xSize)
								break; // Out of bounds.
							else if (j + yy >= 0 && j + yy < ySize)
							{
								// Color in green.
								input[(i + xx)*ySize * 3 + j * 3] = 0; // R
								input[(i + xx)*ySize * 3 + j * 3 + 1] = 255; // G
								input[(i + xx)*ySize * 3 + j * 3 + 2] = 0; // B
							}
						}
					}
				}
			}
		}
	}

	delete[] R;
	delete[] sobelVertical;
	delete[] sobelHorizontal;
	delete[] U_buff;
	delete[] V_buff;
}


