// 水印处理.cpp : 定义控制台应用程序的入口点。
//

#include "Watermark.h"
#include <string.h>
#include <stdio.h>


int main(int argc, char* argv[])
{
	CWatermark wm(44100,20000);
	const int size = 64000*5;
	short* buff = new short[size];
	memset(buff,0,size*2);
	FILE* pf = fopen("test.pcm","rb");
	fread(buff,2,64000*5,pf);
	wm.GenerateWatermark(buff,size);
	fclose(pf);
	pf = fopen("water.pcm","rb");
	fread(buff,2,64000*5,pf);
	bool r =	wm.DetectWatermark(buff,size);
	if (r)
	{
		printf("检测到水印\n");
	}
	else
	{
		printf("没有检测到水印\n");
	}

	fclose(pf);
	return 0;
}
 
