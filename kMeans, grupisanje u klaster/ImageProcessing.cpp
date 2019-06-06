
#include "ImageProcessing.h"
#include "FeatureDetection.h"
#include "ColorSpaces.h"
#include "kMeans.h"
#include <QDebug>

void imageProcessingFun(const QString& progName, QImage* const outImgs, const QImage* const inImgs, const QVector<double>& params)
{
	int X_SIZE = inImgs->width();
	int Y_SIZE = inImgs->height();

	*outImgs = QImage(*inImgs);

	if (progName == "Color kMeans")
	{
		ColorBasedKMeans(outImgs->bits(), outImgs->width(), outImgs->height(), floor(params[0]));
	}
	else if (progName == "Intensity+Position kMeans")
	{
		IntensityPlusPositionBasedKMeans(outImgs->bits(), outImgs->width(), outImgs->height(), floor(params[0]));
	}
	else if (progName == "Color+Position kMeans")
	{
		ColorPlusPositionBasedKMeans(outImgs->bits(), outImgs->width(), outImgs->height(), floor(params[0]));
	}
	else if (progName == "SIFT kMeans")
	{
		SIFTBasedKMeans(outImgs->bits(), outImgs->width(), outImgs->height(), floor(params[0]));
	}
}

