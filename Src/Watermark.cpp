#include "Watermark.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <string.h>

#define TWO_PI          (6.28318530717958647693)
#define FFT_LENGTH		2048
//帧移 ms
#define FFT_SHIFT		10 
#define FFT_WINDOWS		25
#define FFT_AVOID       2

CWatermark::CWatermark(void)
{
	nSamples_ = 44100;
	nFreq_ = 20000;
	nFFT_ = FFT_LENGTH;
	nWindowsShift_ = FFT_SHIFT*nSamples_/1000;
	nWindowsSize_ = FFT_WINDOWS*nSamples_/1000;
}

CWatermark::CWatermark( int nSamples,int nFreq )
{
	nSamples_ = nSamples;
	nFreq_ = nFreq;
	nFFT_ = FFT_LENGTH;
	nWindowsShift_ = FFT_SHIFT*nSamples_/1000;
	nWindowsSize_ = FFT_WINDOWS*nSamples_/1000;
}


CWatermark::~CWatermark(void)
{
}

void CWatermark::GenerateWatermark( short* pdata,int nSize )
{
	static double wt = 0;
	float data = 0;
	float* pOut = new float[nSize];
	float max = 0,ratio = 1.0;
	for (int i = 0 ;i < nSize; i++)
	{
		pOut[i] = pdata[i] + (8000)*sin(wt);		
		wt =  ((TWO_PI * nFreq_)/nSamples_)*i;
		if (wt > TWO_PI)
		{
			wt -= TWO_PI;
		}
		if ( abs(pOut[i]) > max )
		{
			max = abs(pOut[i]);
		}
	}
	// 防止截幅
	if (max > 32000)
	{
		ratio = 32000 / max;
	}
	for (int i = 0 ;i < nSize; i++)
	{
		pdata[i] = pOut[i]*ratio;
	}
	delete [] pOut;
#if 0
	FILE* fp = fopen("water.pcm","wb");
	if (fp)
	{
		fwrite(pdata,2,nSize,fp);
		fclose(fp);
	}
#endif
}

#define SQUARE(a) ((a)*(a))

bool CWatermark::DetectWatermark( const short* pdata,int nSize )
{
	int nWindows = nSize / nWindowsShift_;
	int nlowIndex = nFreq_*nFFT_/nSamples_;
	int nHighIndex = nFreq_*nFFT_/nSamples_+ 10;
	int nNullIndex = nHighIndex + 20 ;
	const int nLowIndex = 50;
	if (nNullIndex >= nFFT_/2 || nWindows < 10 )
	{
		return false;
	}
	float* pfftResult = new float[nWindows];
	float* pfftNull = new float[nWindows];
	float* pfftLow = new float[nWindows];
	for (int i = 0;i<nWindows;i++)
	{
		memset(fftBuf_,0,4096*sizeof(float));
		for (int j=0;j<nWindowsSize_;j++)
		{
			fftBuf_[j] = pdata[i*nWindowsShift_+j]/32678.0;
		}
		preEmphasis_and_hammingWin(fftBuf_,nWindowsSize_);
		doRealfft(fftBuf_,nFFT_,1);
		pfftResult[i] =sqrt( SQUARE(fftBuf_[2*nlowIndex])+SQUARE(fftBuf_[2*nlowIndex+1]));
		pfftResult[i] += sqrt(SQUARE(fftBuf_[2*nHighIndex])+SQUARE(fftBuf_[2*nHighIndex+1]));
		pfftNull[i] = sqrt( SQUARE(fftBuf_[2*nNullIndex])+SQUARE(fftBuf_[2*nNullIndex+1]));
		pfftLow[i] = sqrt( SQUARE(fftBuf_[2*nLowIndex])+SQUARE(fftBuf_[2*nLowIndex+1]));
	}

	float avg0,dev0,avg1,dev1,avg2,dev2;
	// 0序列 均值大 方差小 
	// 1序列 均值小 方差小
	// 2序列 均值中 方差大
	deviation(pfftResult+FFT_AVOID,nWindows-FFT_AVOID*2,avg0,dev0);
	deviation(pfftNull+FFT_AVOID,nWindows-FFT_AVOID*2,avg1,dev1);
	deviation(pfftLow+FFT_AVOID,nWindows-FFT_AVOID*2,avg2,dev2);
	if (avg2 == 0.0){avg2 = 1.0;}
	delete [] pfftResult;
	delete [] pfftNull;
	delete [] pfftLow;
	if ( avg0/avg2 > 0.25 && dev0 < 10 && avg1 < 10 && dev1 < 2)
	{
		return true;
	}
	return false;
}

void CWatermark::preEmphasis_and_hammingWin( float * buff, const int size )
{
	float f = (float)(TWO_PI / (size - 1));
	for (int i = size - 1;i > 0;i--)
	{
		buff[i] -= buff[i - 1] * 0.97;
		buff[i] *= (float)(0.54 - 0.46 * cos(f * i));
	}
	buff[0] *= (float)(1.0 -0.97);
	buff[0] *= (float)(0.54 - 0.46);
}

int CWatermark::dofft (float *data, long nn, int isign)
{
	long n = nn << 1, mmax = 2, m, j = 0, i;
	for (i = 0; i < n - 1; i += 2)
	{
		if (j > i)
		{
			float dum;
			dum = data [j], data [j] = data [i], data [i] = dum;
			dum = data [j+1], data [j+1] = data [i+1], data [i+1] = dum;
		}
		m = n >> 1;
		while (m >= 1 && j + 1 > m)
		{
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	while (n > mmax)
	{
		long istep = 2 * mmax;
		float theta = 2 * M_PI / (isign * mmax);
		float wr, wi, wtemp, wpr, wpi;
		wtemp = sin (0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin (theta);
		wr = 1.0, wi = 0.0;
		for (m = 0; m < mmax - 1; m += 2)
		{
			for (i = m; i < n; i += istep)
			{
				float tempr, tempi;
				j = i + mmax;
				tempr = wr * data [j] - wi * data [j+1], tempi = wr * data [j+1] + wi * data [j];
				data [j] = data [i] - tempr, data [j+1] = data [i+1] - tempi;
				data [i] += tempr, data [i+1] += tempi;
			}
			wtemp = wr, wr = wr * wpr - wi * wpi + wr, wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
	return 0;
}

int CWatermark::doRealfft(float *data, long n, int isign)
{
	long i, i1, i2, i3, i4, np3;
	float c1 = 0.5, c2, h1r, h1i, h2r, h2i;
	float wr, wi, wpr, wpi, wtemp, theta;
	theta = M_PI / (float) (n >> 1);
	if (isign == 1)
	{
		c2 = -0.5;
		dofft(data, n >> 1, 1);
	}
	else
	{
		c2 = 0.5;
		theta = - theta;
	}
	wtemp = sin (0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin (theta);
	wr = 1.0 + wpr;
	wi = wpi;
	np3 = n + 1;
	for (i = 1; i < n >> 2; i++)
	{
		i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i )));
		h1r = c1 * (data [i1] + data [i3]);
		h1i = c1 * (data [i2] - data [i4]);
		h2r = - c2 * (data [i2] + data [i4]);
		h2i = c2 * (data [i1] - data [i3]);
		data [i1] = h1r + wr * h2r - wi * h2i;
		data [i2] = h1i + wr * h2i + wi * h2r;
		data [i3] = h1r - wr * h2r + wi * h2i;
		data [i4] = - h1i + wr * h2i + wi * h2r;
		wr = (wtemp = wr) * wpr - wi * wpi + wr;
		wi = wi * wpr + wtemp * wpi + wi;
	}
	if (isign == 1)
	{
		data [0] = (h1r = data [0]) + data [1];
		data [1] = h1r - data [1];
	}
	else
	{
		data [0] = c1 * ((h1r = data [0]) + data [1]);
		data [1] = c1 * (h1r - data [1]);
		dofft(data, n >> 1, -1);
	}
	return 0;
}

void CWatermark::deviation(float* f, int size,float& average,float& dev )
{
	int i;
	average=0,dev=0;
	for(i=0;i<size;i++)
	{
		average+=f[i]; 
	}

	average/=size;//平均值

	for(i=0;i<size;i++)
	{
		dev+=pow(f[i]-average,2);// 偏离平均数的距离和 
	}

	dev=dev/size;//方差 
	dev=sqrt(dev); //标准差 
}