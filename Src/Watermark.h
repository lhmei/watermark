
#ifndef Watermark_h__
#define Watermark_h__


class CWatermark
{
public:
	CWatermark(void);
	CWatermark(int nSamples,int nFreq);
	~CWatermark(void);

	void GenerateWatermark(short* pdata,int nSize);
	bool DetectWatermark(const short* pdata,int nSize);
protected:
	void preEmphasis_and_hammingWin( float * buff, const int size );
	int dofft (float *data, long nn, int isign);
	int doRealfft(float *data, long n, int isign);
	void deviation(float* f, int size,float& average,float& dev );
private:
	int nSamples_;
	int nFreq_;
	int nWindowsShift_;
	int nWindowsSize_;
	int nFFT_;
	float fftBuf_[4096];
};
#endif // Watermark_h__

