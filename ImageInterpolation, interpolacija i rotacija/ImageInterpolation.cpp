#include "ImageInterpolation.h"
#include <math.h>

//losiji nacin interpolacije
void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
	double fx = ((double)newXSize) / xSize;
	double fy= ((double)newYSize) / ySize;
	int i = 0, j = 0;
	for (i = 0; i < newXSize; i++) {
		for (j = 0; j < newYSize; j++) {
			int newI = round((double)i / fx);
			int newJ = round((double)j / fy);

			if (newI >= xSize) {
				newI = xSize - 1;
			}
			if (newJ >= ySize) {
				newJ = ySize - 1;
			}
			
			output[3*j*newXSize+3*i] = input[3*newJ*xSize+3*newI];
			output[3*j*newXSize + 3 * i+1] = input[3 * newJ*xSize + 3 * newI+1];
			output[3 * j*newXSize + 3 * i + 2] = input[3 * newJ*xSize + 3 * newI + 2];

		}
	}
}
//bolji nacin interpolacije
void bilinearInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
	double fx = ((double)newXSize) / xSize;
	double fy = ((double)newYSize) / ySize;
	int i = 0, j = 0;
	for (i = 0; i < newXSize; i++)
	{
		for (j = 0; j < newYSize; j++)
		{
			int newI = floor((double)i / fx);
			int newJ = floor((double)j / fy);

			if (newI >= xSize-1) {
				newI = xSize - 1;
			}
			if (newJ >= ySize-1) {
				newJ = ySize - 1;
			}

			double b = j / fy - floor(j / fy);
			double a = i / fx - floor(i / fx);
			/*
			y=(1-a)(1-b) X(m,n)+(1-a)b X(m,n)+a(1-b) X(m,n)+abX(m,n)

			*/

			output[3 * i + j * newXSize * 3] =
				(1 - a) * (1 - b) * input[3 * newI + newJ * xSize * 3] +
				(1 - a) * b * input[3 * newI + (newJ + 1) * xSize * 3] +
				a * (1 - b) * input[3 * (newI + 1) + newJ * xSize * 3] +
				a * b * input[3 * (newI + 1) + (newJ + 1) * xSize * 3];
			output[3 * i + 1 + j * newXSize * 3] =
				(1 - a) * (1 - b) * input[3 * newI + 1 + newJ * xSize * 3] +
				(1 - a) * b * input[3 * newI + 1 + (newJ + 1) * xSize * 3] +
				a * (1 - b) * input[3 * (newI + 1) + 1 + newJ * xSize * 3] +
				a * b * input[3 * (newI + 1) + 1 + (newJ + 1) * xSize * 3];
			output[3 * i + 2 + j * newXSize * 3] =
				(1 - a) * (1 - b) * input[3 * newI + 2 + newJ * xSize * 3] +
				(1 - a) * b * input[3 * newI + 2 + (newJ + 1) * xSize * 3] +
				a * (1 - b) * input[3 * (newI + 1) + 2 + newJ * xSize * 3] +
				a * b * input[3 * (newI + 1) + 2 + (newJ + 1) * xSize * 3];
		}
	}
}
	

//rotacija slike
void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
	
	 for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            double di = i - m;
            double dj = j - n;
         
            double theta = angle *  3.14159265359 / 180;

            int newI = round(i* cos(theta) - j * sin(theta) -m*cos(theta) + n*sin(theta) + m);
            int newJ = round(j* cos(theta) + i* sin(theta) - m*sin(theta) -n* cos(theta) + n);

            if (newI >= 0 && newI < xSize && newJ >= 0 && newJ < ySize)
            {
                output[3 * i + j * xSize * 3] = input[3 * newI + newJ * xSize * 3];
                output[3 * i + 1 + j * xSize * 3] = input[3 * newI + 1 + newJ * xSize * 3];
                output[3 * i + 2 + j * xSize * 3] = input[3 * newI + 2 + newJ * xSize * 3];
            }
            else
            {
                output[3 * i + j * xSize * 3] = 0;
                output[3 * i + 1 + j * xSize * 3] = 0;
                output[3 * i + 2 + j * xSize * 3] = 0;
            }
        }
    }
	

}

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */

	for (int i = 0; i < xSize; i++)
	{
		for (int j = 0; j < ySize; j++)
		{
			double di = i - m;
			double dj = j - n;

			double theta = angle *  3.14159265359 / 180;

			double rotI = i* cos(theta) - j * sin(theta) - m*cos(theta) + n*sin(theta) + m;
			double rotJ = j* cos(theta) + i* sin(theta) - m*sin(theta) - n* cos(theta) + n;

			int newI = floor(rotI);
			int newJ = floor(rotJ);

			if (newI >= xSize - 1) {
				newI = xSize - 1;
			}
			if (newJ >= ySize - 1) {
				newJ = ySize - 1;
			}

			double a = rotJ - newJ;
			double b = rotI - newI;

			if (newI >= 0 && newI + 1 < xSize && newJ >= 0 && newJ + 1 < ySize)
			{
				output[3 * i + j * xSize * 3] =
					(1 - a) * (1 - b) * input[3 * newI + newJ * xSize * 3] +
					(1 - a) * b * input[3 * (newI + 1) + newJ * xSize * 3] +
					a * (1 - b) * input[3 * newI + (newJ + 1) * xSize * 3] +
					a * b * input[3 * (newI + 1) + (newJ + 1) * xSize * 3];
				output[3 * i + 1 + j * xSize * 3] =
					(1 - a) * (1 - b) * input[3 * newI + 1 + newJ * xSize * 3] +
					(1 - a) * b * input[3 * (newI + 1) + 1 + newJ * xSize * 3] +
					a * (1 - b) * input[3 * newI + 1 + (newJ + 1) * xSize * 3] +
					a * b * input[3 * (newI + 1) + 1 + (newJ + 1) * xSize * 3];
				output[3 * i + 2 + j * xSize * 3] =
					(1 - a) * (1 - b) * input[3 * newI + 2 + newJ * xSize * 3] +
					(1 - a) * b * input[3 * (newI + 1) + 2 + newJ * xSize * 3] +
					a * (1 - b) * input[3 * newI + 2 + (newJ + 1) * xSize * 3] +
					a * b * input[3 * (newI + 1) + 2 + (newJ + 1) * xSize * 3];
			}
			else
			{
				output[3 * i + j * xSize * 3] = 0;
				output[3 * i + 1 + j * xSize * 3] = 0;
				output[3 * i + 2 + j * xSize * 3] = 0;
			}
		}
	}
}


