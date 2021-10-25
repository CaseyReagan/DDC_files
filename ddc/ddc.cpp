// ddc.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <string>
#include <iostream>  
#include <math.h>
//#include <ipp.h>
#include <ipps.h>
#include <cmath>
//#include <ddc.h>

//using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	int i=0,j=0,k=0,h=0;

	const int l = 129;				//length
	const int pl = 16384;			//pipe_length
	const double pi = 3.14159265358979323846264338328;
	const double E = 2.71828182845904523536029;

	double s_freq = 186.6;			//source_frequency 186.6MHz
	double t_freq1 = 46.69;			//target_frequency1 46.69MHz
	double t_samp_rate = 8.0;		//taiget sampling rate 8MHz
	int res = 1,count=0;
	int insert_mult = pow(2,floor(log(s_freq / t_samp_rate) / log(2)));    //插值倍数 2 的取整 log(x)/log(2) (in C++) = log2(x) (in math)次幂

	//double fw1 = 16/186.6;		//used in insert data 
	//double fw2 = 1/8.0;
	double fw1 = insert_mult/s_freq;		//used in insert data 
	double fw2 = 1/t_samp_rate;
	double tl = fw1;
	double d_w = fw2 - fw1;			//宽度差
	double u = 0, hm=0, u1=0;
	double im[10] = {0};
	bool jump = 0;

	char *addr = "F:\\program_files\\ddc\\sig_source.dat";
	char *addr2 = "F:\\program_files\\ddc\\rcosfir.txt";

	__int16 source[pl+l];					//source signal data
//	__int16 *source = new __int16[pl+l];	//source signal data
	Ipp64f *ip11 = new Ipp64f[pl+l];		//变频I路
	Ipp64f *qp11 = new Ipp64f[pl+l];		//变频q路
	Ipp64f *ip21 = new Ipp64f[pl+l];		//变频I路
	Ipp64f *qp21 = new Ipp64f[pl+l];		//变频q路
	Ipp64f *ip22 = new Ipp64f[pl+l];		//变频I路
	Ipp64f *qp22 = new Ipp64f[pl+l];		//变频q路
	Ipp64f *ip1 = new Ipp64f[2*pl+l];					//变频I路综合观测1
//	Ipp64f *qp1	= new Ipp64f[2*pl+l];					//变频q路综合观测1
	Ipp64f *ip31 = new Ipp64f[pl+l];		//CIC I路
	Ipp64f *qp31 = new Ipp64f[pl+l];		//CIC q路
	Ipp64f *ip32 = new Ipp64f[pl/16];		//CIC I路抽值
	Ipp64f *qp32 = new Ipp64f[pl/16];		//CIC q路抽值
	Ipp64f *ip33 = new Ipp64f[pl/16];		//CIC I路
	Ipp64f *qp33 = new Ipp64f[pl/16];		//CIC q路
	Ipp64f *ip41 = new Ipp64f[pl/16];		//CIC q路

	Ipp64f y1=0,y2=0,y3=0;					//CIC参数
	Ipp64f x1=0,x2=0,x3=0;					//CIC参数
	Ipp64f y1m=0,y2m=0,y3m=0;
	Ipp64f x1m=0,x2m=0,x3m=0;

	Ipp64f y4=0,y5=0,y6=0;					//CIC参数
	Ipp64f x4=0,x5=0,x6=0;					//CIC参数

	Ipp64f ym1=0,ym2=0,ym3=0;				//CIC参数
	Ipp64f xm1=0,xm2=0,xm3=0;				//CIC参数

	Ipp64f monitor;

	Ipp64f rFreq = 1 / s_freq;				 //IPP滤波所需参数
	int tapslen = l-1;
	int numIters = pl + l;
	const Ipp64f *ipp11 = ip11;				 //pointer point to the array ip11
	const Ipp64f *qpp11 = qp11;				 //pointer point to the array qp11
	IppsFIRState_64f *pState;
	Ipp64f *FIRDst = ippsMalloc_64f((pl + l)*sizeof(Ipp64f));		// remember to ippsFree()
	Ipp64f *taps = ippsMalloc_64f(tapslen * sizeof(Ipp64f));
	Ipp64f *pDL = ippsMalloc_64f(tapslen * sizeof(Ipp64f));
	ippsZero_64f(pDL,tapslen);
	/****  computes tapsLen coefficients for lowpass FIR filter  ****/
	ippsFIRGenLowpass_64f(rFreq, taps, tapslen, ippWinHamming, ippTrue);
	ippsFIRInitAlloc_64f(&pState, taps,tapslen, pDL);

	int tapslen2 = 48;
	IppsFIRState_64f *pState2;				//rise cos fir
	Ipp64f *taps2 = ippsMalloc_64f(tapslen2 * sizeof(Ipp64f));
	Ipp64f *pDL2 = ippsMalloc_64f(tapslen2 * sizeof(Ipp64f));
	ippsZero_64f(pDL2,tapslen2);

	FILE *fp2;
	fp2 = fopen(addr2,"rb");
	if(fp2 == NULL)
		printf("file open error.");
	fread(taps2,sizeof(double),tapslen2,fp2);
	fclose(fp2);

	ippsFIRInitAlloc_64f(&pState2, taps2,tapslen2, pDL2);

	FILE *fp;
	fp = fopen(addr,"rb");
	if(fp == NULL)
		printf("file open error.");

	/**** get the source signal data ****/
//	fread(source,sizeof(__int16),pl+l,fp);
	while(res)
	{
		if(fread(source,sizeof(__int16),pl+l,fp) == pl+l)             //读取数据的两段完全正确
		{
			res = 1;

			/**** 变频操作 ****/
			for(i=j,k=0;i<j+pl+l,k<pl+l;i++,k++)
			{
				ip11[k] = static_cast<Ipp64f> (source[k] * cos(2 * pi * (t_freq1/s_freq) * (i+1)));
				qp11[k] = static_cast<Ipp64f> (source[k] * sin(2 * pi * (t_freq1/s_freq) * (i+1)));
			}
			j =  i - l;

/*			for(i=count*pl,k=0;i<(count+1)*pl+l,k<pl+l;i++,k++)
			{
				ip11[k] = static_cast<Ipp64f> (source[k] * cos(2 * pi * (t_freq1/s_freq) * (i+1)));
				qp11[k] = static_cast<Ipp64f> (source[k] * sin(2 * pi * (t_freq1/s_freq) * (i+1)));
			}
*/

			//验证变频结果是否正确
			//结果正确
			for(i=count*pl,k=0;k<pl;i++,k++)
			{
				ip1[i] = ip11[k];
			}


			/**** 滤波 ****/
			//filter an input vector+*-------------------
			ippsFIR_64f(ipp11, FIRDst, numIters, pState);

			for(k=0;k<pl+l;k++)
			{
				ip21[k] = FIRDst[k];
			}
			ippsFIR_64f(qpp11, FIRDst, numIters, pState);
			for(k=0;k<pl+l;k++)
			{
				qp21[k] = FIRDst[k];
			}

			/****  阶段综合结果检测，送给matlab  ****/
			/*******************  滤波组成的数据结果经验证完全正确  *******************/
/*			if(count<1)
			{
				for(k=0;k<pl+l;k++)
				{
					ip1[k] = ip21[k];
				}
			}
			else
			{
				for(i=pl+l+(count-1)*pl,k=l;k<pl+l;i++,k++)
				{
					ip1[i] = ip21[k];
				}
			}
*/
			/**** 用上一轮最后l个数替换本轮数据头l个数 ****/
			if(count > 0)						//第一轮的结果不用改变,数据头有一点不准确可以忽略，将第一段结果存在ip22里
			{
				for(k=0;k<l;k++)				//第二轮开始，将上一轮数据最后l位，替换新一轮数据的前l位
				{
					ip21[k] = ip22[k+pl];
					qp21[k] = qp22[k+pl];
				}
			}

			for(k=0;k<pl+l;k++)					//将本轮结果保存在ip22中，用来给下一轮用
			{
				ip22[k] = ip21[k];
				qp22[k] = qp22[k];
			}

			/*******************  滤波组成的数据结果至这一步,经验证完全正确  *******************/
/*			for(i=count*pl,k=0;k<pl;i++,k++)
			{
				ip1[i] = ip21[k];			  //"filter_check3.dat"
			}
*/

			/****  CIC  ****/
			/****  CIC first step  ****/
			y1 = y1m;
			y2 = y2m;
			y3 = y3m;

			x1 = x1m;
			x2 = x2m;
			x3 = x3m;

			for(i=0;i<pl+l;i++)
			{
				y1 = ip21[i] + y1;
				y2 = y1 + y2;
				y3 = y2 + y3;
				ip31[i] = y3;

				x1 = qp21[i] + x1;
				x2 = x1 + x2;
				x3 = x2 + x3;
				qp31[i] = x3;

				if(i == pl-1)
				{
					y1m = y1;
					y2m = y2;
					y3m = y3;

					x1m = x1;
					x2m = x2;
					x3m = x3;
				}
			}

/*			for(i=count*pl,k=0;k<pl;i++,k++)
			{
				ip1[i] = ip31[k];				//"filter_check4.dat"
			}
*/

			/****  CIC second step 抽值 ****/
			for(i=0,k=0;i<pl/insert_mult;i++,k+=16)
			{
				ip32[i] = ip31[k];
				qp32[i] = qp31[k];
			}

			/******************* CIC组成的数据结果至这一步,经验证完全正确  *******************/
/*			for(k=count*pl/samp_mult,i=0;k<(count+1)*pl/samp_mult;i++,k++)
			{
				ip1[k] = ip32[i];              //"filter_check5.dat"
			}
*/

			/****  CIC final step  ****/
			/******************* CIC组成的数据结果至这一步,经验证完全正确  *******************/
			for(i=0;i<pl/insert_mult;i++)
			{
				ym1 = ip32[i] - y4;
				y4 = ip32[i];
				ym2 = ym1 - y5;
				y5 = ym1;
				ym3 = ym2 - y6;
				y6 = ym2;
				ip33[i] = ym3;
				
				xm1 = qp32[i] - x4;
				x4 = qp32[i];
				xm2 = xm1 - x5;
				x5 = xm1;
				xm3 = xm2 - x6;
				x6 = xm2;
				qp33[i] = xm3;
			}

/*			for(k=count*pl/samp_mult,i=0;k<(count+1)*pl/samp_mult;i++,k++)
			{
				ip1[k] = ip41[i];				//"filter_check7.dat"
			}
*/
			/****  insert data  ****/
			if(count == 0)										//第一段特殊处理
			{
				for(i=0,h=1;h<pl/insert_mult;h++)
				{
					if(h > 1)
					{
						u = u + d_w;
						if(u > tl)
						{
							u = u - tl;
							h = h + 1;
						}
					}

					if(h < pl/insert_mult)
					{
						ip41[i] = u / tl * (ip33[h]-ip33[h-1]) + ip33[h-1];		//"filter_check7.dat"
						i = i + 1;
					}

					if(h == pl/insert_mult)
					{
						jump = 1;
					}

					if(h >=  (pl/insert_mult)-1)
					{
						hm = ip33[pl/insert_mult -1];				//记录最后一点的值
						u1 = u;									//记录u值
						im[1] = i;
					}
				}
			}
			else if(count > 0)									//第二段数据开始处理
			{
				u = u1;
				for(i=0,h=0;h<pl/insert_mult;h++)
				{
					if(jump)
					{
						ip41[0] = u / tl * (ip33[0]-hm) + hm;		//"filter_check7.dat"
						jump = 0;

						i++;
					}
					else if(h == 0)
					{
						u = u + d_w;
						if(u > tl)
						{
							u = u - d_w;
						}
						else
						{
							ip41[0] = u / tl * (ip33[0]-hm) + hm;

							i++;
						}
					 }

					if(h > 0 && h < pl/insert_mult)
					{
						u = u + d_w;
						if(u > tl)
						{
							u = u - tl;
							h = h + 1;
						}

						ip41[i] = u / tl * (ip33[h]-ip33[h-1]) + ip33[h-1];		//"filter_check7.dat"
						i = i + 1;
					}

					if(h == pl/insert_mult)
					{
						jump = 1;
					}

					if(h >=  (pl/insert_mult)-1)
					{
						hm = ip33[pl/insert_mult -1];				//记录最后一点的值
						u1 = u;									//记录u值
						im[count+1] = im[count] + i;
					}
				}
			}

			//for(k=im[count],i=0;k<im[count+1];i++,k++)
			//{
			//	ip1[k] = ip41[i];				//"filter_check7.dat"
			//}

			count++;			//计数

			fseek(fp,-(2*l),SEEK_CUR);
		}
		else
		{
			res = 0;
			break;
		}
	}

	//matlab check
	//FILE *fid2=fopen("filter_check7.dat","wb");
	//fwrite(ip1,sizeof(Ipp64f),im[2],fid2);
	//fclose(fid2);

	FILE *fid2=fopen("filter_check1.dat","wb");
	fwrite(ip1,sizeof(Ipp64f),2*pl,fid2);
	fclose(fid2);

	fclose(fp);

	ippsFree(FIRDst);
	ippsFree(taps);
	ippsFree(pDL);

	ippsFree(taps2);
	ippsFree(pDL2);

	//delete []source;
	//source = 0;

	delete []ip11;
	ip11 = 0;

	delete []qp11;
	qp11 = 0;

	delete []ip21;
	ip21 = 0;

	delete []qp21;
	qp21 = 0;

	delete []ip1;
	ip1 = 0;

//	delete []qp1;
//	qp1 = 0;

	delete []ip22;
	ip22 = 0;

	delete []qp22;
	qp22 = 0;

	delete []ip31;
	ip31 = 0;

	delete []qp31;
	qp31 = 0;

	delete []ip32;
	ip32 = 0;

	delete []qp32;
	qp32 = 0;

	delete []ip33;
	ip33 = 0;

	delete []qp33;
	qp33 = 0;

	delete []ip41;
	ip41 = 0;

	system("pause");
	return 0;
}


/*   
	//useless
int get_opt_times(void)
{
	char *addr = "F:\\program_files\\ddc\\sig_source.dat";
//	ifstream file(addr);
	
	return 0;
}
*/
