#include "ImageFilter.h"


void convolve2D(uchar Y_buff[], int xSize, int ySize, double filterCoeff[], int N)
{
	//TO DO
	int newXSize = xSize + (N - 1); int newYSize = ySize + (N - 1);
	int delta = (N - 1) / 2;
	uchar *Y_ext = new uchar[(xSize + 2 * delta)*(ySize + 2 * delta)];
	extendBorders(Y_buff, xSize, ySize, Y_ext, delta);

	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {
			double accum = 0;
			for (int m = 0; m < N; m++)
			{
				for (int n = 0; n < N; n++)
				{
					accum += Y_ext[(j + n) * newXSize + i + m] * filterCoeff[(N - n) * N - m - 1];
				}
			}
			if (accum > 255.0) Y_buff[j * xSize + i] = 255;
			else if (accum < 0.0) Y_buff[j * xSize + i] = 0; else Y_buff[j * xSize + i] = floor(accum + 0.5);
		}
	}
}

void extendBorders(uchar input[], int xSize, int ySize, uchar output[], int delta)
{
	//TO DO
	for (int i = 0; i < delta; i++)
	{
		memcpy(&output[i * (xSize + delta * 2) + delta], &input[0], xSize);
		memcpy(&output[(ySize + delta + i)*(xSize + delta * 2) + delta], &input[(ySize - 1) * xSize], xSize);
	}

	for (int i = 0; i < (ySize + delta * 2); i++)
	{
		memset(&output[i * (xSize + delta * 2)], output[i * (xSize + delta * 2) + delta], delta);
		memset(&output[i * (xSize + delta * 2) + delta + xSize], output[i * (xSize + delta * 2) + xSize + delta - 1], delta);
	}
	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {
			output[(i + delta) + (xSize + 2 * delta)*(j + delta)] = input[i + xSize*j];
		}
	}

}

void performNFFilter(uchar input[], int xSize, int ySize)
{
	//TO DO
	// 1 1 1 
	// 1 1 1    *1/9
	// 1 1 1
	double filterCoefs[9] = { 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9 };
	convolve2D(input, xSize, ySize, filterCoefs, 3);
}

void performVFFilter(uchar input[], int xSize, int ySize)
{
	//TO DO
	//  0 -1  0
	// -1  4 -1
	//  0 -1  0
	double filterCoefs[9] = { 0,-1,0,-1,4,-1,0,-1,0 };
	convolve2D(input, xSize, ySize, filterCoefs, 3);
}

void performSuccessiveNFFilter(uchar input[], int xSize, int ySize, int stages)
{
	//TO DO
	double filterCoefs[9] = { 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9 };
	for (int i = 0; i < stages; i++) {
		convolve2D(input, xSize, ySize, filterCoefs, 3);
	}
}
//za detekciju ivica
void performSobelEdgeDetection(uchar input[], int xSize, int ySize, uchar threshold)
{
	//TO DO
	int N = 3;
	int delta = (N - 1) / 2;
	int newXSize = xSize + (N - 1);
	int newYSize = ySize + (N - 1);
	uchar *Y_ext = new uchar[(xSize + 2 * delta) *(ySize + 2 * delta)];
	extendBorders(input, xSize, ySize, Y_ext, delta);


	double verticalCoefs[9] = { -1.0 / 4, 0, 1.0 / 4, -2.0 / 4, 0, 2.0 / 4, -1.0 / 4, 0, 1.0 / 4 };
	double horizontalCoefs[9] = { -1.0 / 4, -2.0 / 4, -1.0 / 4, 0, 0, 0, 1.0 / 4, 2.0 / 4, 1.0 / 4 };


	for (int i = 0; i < xSize; i++)
	{
		for (int j = 0; j < ySize; j++)
		{
			double hAccum = 0;
			double vAccum = 0;
			for (int m = 0; m < N; m++) {
				for (int n = 0; n < N; n++) {
					hAccum += Y_ext[(j + n) * newXSize + i + m] * horizontalCoefs[(N - n) * N - m - 1];
					vAccum += Y_ext[(j + n) * newXSize + i + m] * verticalCoefs[(N - n) * N - m - 1];

				}
			}
			if (sqrtf(pow(hAccum,2)+pow(vAccum,2)) >= threshold)
			{
				input[i + xSize * j] = 255;
			}
			else
			{
				input[i + xSize * j] = 0;
			}
		}
	}
}

void performNFplusSobelEdgeDetection(uchar input[], int xSize, int ySize, int stages, uchar threshold)
{
	//TO DO
	performSuccessiveNFFilter(input, xSize, ySize, stages);
	performSobelEdgeDetection(input, xSize, ySize, threshold);
}
