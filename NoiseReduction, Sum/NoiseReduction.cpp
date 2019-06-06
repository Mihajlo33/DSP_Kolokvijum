#include "NoiseReduction.h"
#include "ImageFilter.h"

#include <vector>
#include <algorithm>

//popravljanje suma, dosta lose
void performMovingAverage (uchar input[], int xSize, int ySize, int N)
{
	//TO DO
	//kako smo imali 3x3  =9  -->  1/9* [ 1 1 1, 1 1 1,  1 1 1]
	//5x5=25-->1/25*[1 1 1 1 1,1 1 1 1 1,1 1 1 1 1,1 1 1 1 1,1 1 1 1 1]

	double *filter = new  double[N*N];
	for (int i = 0; i < N*N; i++) {
		filter[i] = 1.0 / N*N;
	}


	convolve2D(input, xSize, ySize, filter, N);

	
}

void calculateGaussKernel(double kernel[], int N, double sigma)
{
	double C = 0;
	for(int n = 0; n < N; n++)
    {
        for(int k = 0; k < N; k++)
        {
            kernel[n*N+k] = exp( -((n - N/2)*(n - N/2) + (k - N/2)*(k - N/2)) / (2 * sigma * sigma));
			C += kernel[n*N+k];
		}
	}

	C = 1.0 / C;

	for(int n = 0; n < N; n++)
    {
        for(int k = 0; k < N; k++)
        {
            kernel[n*N+k] *= C;
		}
	}
}
//sredjuje Gausov sum, ali ne impulsni
void performGaussFilter (uchar input[], int xSize, int ySize, int N, double sigma)
{
	//TO DO
	double *kernel = new double[N*N];
	//{ 0 };  //prazna
	calculateGaussKernel(kernel, N, sigma);
	convolve2D(input,xSize, ySize, kernel, N);



}
//sredjuje impulsni sum
void performMedianFilter(uchar input[], int xSize, int ySize, int N)
{
	//TO DO
	//Voditi raèuna da je neophodno izvršiti proširenje slike pre blokovske obrade.
	//n/2+1

	uchar* extended = new uchar[(xSize + N - 1) * (ySize + N- 1)];
	extendBorders(input, xSize, ySize, extended, N / 2);
	//N=3;
	//N=5;
	//N = 7;
	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {
			int* median = new int[N * N];

			for (int k = 0; k <N; k++) {
				for (int m = 0; m < N; m++) {
					median[k + N * m] = extended[(i + k) + (j + m) * (xSize + N - 1)];
				}
			}

			for (int k = 0; k <N * N - 1; k++) {
				for (int m = k; m < N * N; m++) {
					if (median[k] > median[m]) {
						int t = median[k];
						median[k] = median[m];
						median[m] = t;
					}
				}
			}
			input[i + j * xSize] = median[N * N / 2];
			delete median;
		}
	}

	delete extended;
}
