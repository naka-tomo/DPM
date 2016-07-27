#pragma once
/*******************************************************************

GaussianTable.h.h		developed by naka_t	2011.02.23

	CRPでのテーブルクラス

  Copyright (C) 2011  naka_t <naka_t@apple.ee.uec.ac.jp>
 *******************************************************************/
#include <vector>

//OpenCVヘッダ
#include <cv.h>
#include <cxcore.h> 
#include <highgui.h>

//OpenCVライブラリ
#pragma comment (lib, "cv210.lib") 
#pragma comment( lib, "cxcore210.lib")
#pragma comment( lib, "highgui210.lib")


// ガウス分布のパラメタ構造体
struct GaussianParams
{
	GaussianParams(int dim);
	GaussianParams(){ dim = 0; };
	
	void Create( int dim_ );				// メモリ確保
	void Init( double **data , int num );	// 与えられたデータから初期値を決定
	GaussianParams Clone();					// 新たなメモリ領域を確保しコピー

	// ガウス分布のパラメタ
	cv::Mat X;
	cv::Mat C;
	cv::Mat m;
	cv::Mat S;
	int r;
	int nu;
	int N;
	int dim;
};


class CGaussianTable
{
public:
	CGaussianTable(GaussianParams &init);
	~CGaussianTable(void);

	void AddData( double *data );				// テーブルにデータを追加
	void DeleteData( double *data );			// テーブルからデータを削除
	int GetNumData(){ return m_param.N; }		// テーブルに座っている人の人数を取得
	double CalcLogLikilihood( double *data );	// データがこのテーブルに座る対数尤度を計算

protected:
	GaussianParams m_param;		// テーブルのパラメタ
	GaussianParams m_init;		// 初期値

	double CalcLogPx( double *addData=NULL);
	void AddData( GaussianParams &p , double *data );
};
