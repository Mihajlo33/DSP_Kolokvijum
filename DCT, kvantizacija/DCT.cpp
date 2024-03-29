#include "DCT.h"
#include "NxNDCT.h"
#include "stdio.h"

/* perform DCT */ /*napravi T matricu, prosiri sliku, uzme svaki blok i nad njim radi DCT i iscrtava ih*/
void performDCT(uchar Y_buff[], int xSize, int ySize, int N)
{
	/* Create NxN buffer for one input block */
	uchar* inBlock = new uchar[N*N];

	/* Generate DCT kernel */
	double* DCTKernel = new double[N*N];
	GenerateDCTmatrix(DCTKernel, N);

	/* Create buffer for DCT coefficients */
	short* dctCoeffs = new short[N*N];

	/* Extend image borders */
	int xSize2, ySize2;
	uchar* Y_buff2;
	extendBorders(Y_buff, xSize , ySize, N, &Y_buff2, &xSize2, &ySize2);

	for (int y = 0; y < ySize2/N; y++)
	{
		for (int x = 0; x < xSize2/N; x++)
		{	

			/* Fill input block buffer */
			for (int j=0; j<N; j++)
			{
				for (int i=0; i<N; i++)
				{
					inBlock[j*N+i] = Y_buff2[(N*y+j)*(xSize2)+N*x+i];
				}
			}

			/* Invoke DCT */
			DCT(inBlock, dctCoeffs, N, DCTKernel);	

			/* Prepere DCT coefficients for drawing and put them in Y_buff */
			for(int i = 0; i < N*N; i++)
			{
				dctCoeffs[i] = abs(dctCoeffs[i]);
				if(dctCoeffs[i] < 0)
				{
					dctCoeffs[i] = 0;
				}
				else if(dctCoeffs[i] > 255)
				{
					dctCoeffs[i] = 255;
				}

				inBlock[i] = dctCoeffs[i];
			}

			/* Write output values to output image */
			for (int j=0; j<N; j++)
			{
				for (int i=0; i<N; i++)
				{
					Y_buff2[(N*y+j)*(xSize2)+N*x+i] = inBlock[j*N+i];
				}
			}
		}
	}

	cropImage(Y_buff2, xSize2, ySize2, Y_buff, xSize, ySize);

	/* Delete used memory buffers coefficients */
	delete[] dctCoeffs;
	delete[] inBlock;
	delete[] DCTKernel;
	delete[] Y_buff2;
}

/* perform DCT and IDCT */
void performDCTandIDCT(uchar* Y_buff, int xSize, int ySize, int N)
{
	/* Create NxN buffer for one input block */
	uchar* inBlock = new uchar[N*N];

	/* Generate DCT kernel */
	double* DCTKernel = new double[N*N];
	GenerateDCTmatrix(DCTKernel, N);

	/* Create buffer for DCT coefficients */
	short* dctCoeffs = new short[N*N];

	/* Extend image borders */
	int xSize2, ySize2;
	uchar* Y_buff2;
	extendBorders(Y_buff, xSize , ySize, N, &Y_buff2, &xSize2, &ySize2);

	for (int y = 0; y < ySize2/N; y++)
	{
		for (int x = 0; x < xSize2/N; x++)
		{	

			/* Fill input block buffer */
			for (int j=0; j<N; j++)
			{
				for (int i=0; i<N; i++)
				{
					inBlock[j*N+i] = Y_buff2[(N*y+j)*(xSize2)+N*x+i];
				}
			}

			/* Invoke DCT */
			DCT(inBlock, dctCoeffs, N, DCTKernel);	



			/* Invoke IDCT */
			IDCT(dctCoeffs, inBlock, N, DCTKernel);	

			/* Write output values to output image */
			for (int j=0; j<N; j++)
			{
				for (int i=0; i<N; i++)
				{
					Y_buff2[(N*y+j)*(xSize2)+N*x+i] = inBlock[j*N+i];
				}
			}
		}
	}

	cropImage(Y_buff2, xSize2, ySize2, Y_buff, xSize, ySize);

	/* Delete used memory buffers coefficients */
	delete[] dctCoeffs;
	delete[] inBlock;
	delete[] DCTKernel;
	delete[] Y_buff2;
}

/* perform masking of DCT coeffs */
void performMaskDCTCoeffs(uchar input[], int xSize, int ySize, int N, int B)
{
	//TO DO
	/* Create NxN buffer for one input block */
	uchar* inBlock = new uchar[N*N];

	/* Generate DCT kernel */
	double* DCTKernel = new double[N*N];
	GenerateDCTmatrix(DCTKernel, N);
	/* Create buffer for DCT coefficients */
	short* dctCoeffs = new short[N*N];

	/* Extend image borders */
	int xSize2, ySize2;
	uchar* Y_buff2;
	//extendBorders(Y_buff, xSize, ySize, N, &Y_buff2, &xSize2, &ySize2);
	extendBorders(input, xSize, ySize, N, &Y_buff2, &xSize2, &ySize2);

	for (int y = 0; y < ySize2 / N; y++)
	{
		for (int x = 0; x < xSize2 / N; x++)
		{

			/* Fill input block buffer */
			for (int j = 0; j<N; j++)
			{
				for (int i = 0; i<N; i++)
				{
					inBlock[j*N + i] = Y_buff2[(N*y + j)*(xSize2)+N*x + i];
				}
			}

			/* Invoke DCT */
			DCT(inBlock, dctCoeffs, N, DCTKernel);

		//////////////////////////////////////
		//maska za gornjih levih BxB bita, jer nose najvecu informaciju o slici, pokazuje da su oni najbitniji
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					if (i > B || j > B) { //veca od velicine bloka
						//output[i*N+j] = floor(sum+0.5);
						dctCoeffs[i*N + j] = 0;
					}
				}
			}

		/////////////////////////////////////
		/* Invoke IDCT */
			IDCT(dctCoeffs, inBlock, N, DCTKernel);

			/* Write output values to output image */
			for (int j = 0; j<N; j++)
			{
				for (int i = 0; i<N; i++)
				{
					Y_buff2[(N*y + j)*(xSize2)+N*x + i] = inBlock[j*N + i];
				}
			}
		}
	}
	
	/////////////////////////////////////////////////
	//izvr�iti maskiranje podmatrice BxB u gornjem levom uglu, tako da sve ostale vrednosti u bloku budu 0 
	double MSE = 0;
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < N; i++)
		{
			MSE += pow(input[i*xSize + j] - Y_buff2[i*xSize2 + j], 2);
		}
	}
	MSE /= xSize*ySize;

	/////////////////////////////////////////////////
	/* Delete used memory buffers coefficients */
	cropImage(Y_buff2, xSize2, ySize2, input, xSize, ySize);

	delete[] dctCoeffs;
	delete[] inBlock;
	delete[] DCTKernel;
	delete[] Y_buff2;
}



/* 8x8 Quantization matrix */
uchar quantizationMatrix[64] = 
{
    16, 11, 10, 16, 24, 40, 51, 61,
    12, 12, 14, 19, 26, 58, 60, 55,
    14, 13, 16, 24, 40, 57, 69, 56,
    14, 17, 22, 29, 51, 87, 80, 62,
    18, 22, 37, 56, 68, 109, 103, 77,
    24, 35, 55, 64, 81, 104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99
};

/* perform quantization of DCT coeffs using Q matrix*/
void performDCTQuantization(uchar input[], int xSize, int ySize)
{
	//TO DO
	int N = 8;
	/* Create NxN buffer for one input block */
	uchar* inBlock = new uchar[N*N];

	/* Generate DCT kernel */
	double* DCTKernel = new double[N*N];
	GenerateDCTmatrix(DCTKernel, N);

	/* Create buffer for DCT coefficients */
	short* dctCoeffs = new short[N*N];

	/* Extend image borders */
	int xSize2, ySize2;
	uchar* Y_buff2;
	extendBorders(input, xSize, ySize, N, &Y_buff2, &xSize2, &ySize2);

	for (int y = 0; y < ySize2 / N; y++)
	{
		for (int x = 0; x < xSize2 / N; x++)
		{

			/* Fill input block buffer */
			for (int j = 0; j<N; j++)
			{
				for (int i = 0; i<N; i++)
				{
					inBlock[j*N + i] = Y_buff2[(N*y + j)*(xSize2)+N*x + i];
				}
			}

			/* Invoke DCT */
			DCT(inBlock, dctCoeffs, N, DCTKernel);
			//////////////////////////////////////////////////////////////////////////////
			for (int j = 0; j < N; j++)
			{
				for (int i = 0; i < N; i++)
				{
					dctCoeffs[i*N + j] = floor(dctCoeffs[i*N + j] / quantizationMatrix[i*N + j]);
					dctCoeffs[i*N + j] *= quantizationMatrix[i*N + j];
				}
			}



			//////////////////////////////////////////////////////////////////////////////
			/* Invoke IDCT */
			IDCT(dctCoeffs, inBlock, N, DCTKernel);

			/* Write output values to output image */
			for (int j = 0; j<N; j++)
			{
				for (int i = 0; i<N; i++)
				{
					Y_buff2[(N*y + j)*(xSize2)+N*x + i] = inBlock[j*N + i];
				}
			}
		}
	}

	
	/////////////////////////////////////////////////
	//greska prilikom kvantizacije
												
	double MSE = 0;
	for (int j = 0; j < xSize; j++)
	{
			for (int i = 0; i < ySize; i++)
		{
			MSE += pow(input[i*xSize + j] - Y_buff2[i*xSize2 + j], 2);
		}
			
	}
	
	MSE /= xSize*ySize;
	printf("MSE= %2.f", MSE);
												
												


	/////////////////////////////////////////////////
	cropImage(Y_buff2, xSize2, ySize2, input, xSize, ySize);

	/* Delete used memory buffers coefficients */
	delete[] dctCoeffs;
	delete[] inBlock;
	delete[] DCTKernel;
	delete[] Y_buff2;

}