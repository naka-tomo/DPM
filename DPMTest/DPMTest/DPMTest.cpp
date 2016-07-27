
/*******************************************************************

DPMTest.cpp		developed by naka_t	2011.02.01

	DPMクラスの使用例

  Copyright (C) 2011  naka_t <naka_t@apple.ee.uec.ac.jp>
 *******************************************************************/
#include "stdafx.h"
#include "MultinomialDPM.h"
#include "utility.h"

int _tmain(int argc, _TCHAR* argv[])
{
	CMultinomialDPM dpm;	// DPMクラス
	int num, dim;			// データ数と次元

	// データの読み込み
	double **data = LoadMatrix<double>( dim , num , "sample.txt" );

	// データを渡す
	dpm.SetData( data , num , dim );

	// パラメタの更新
	for(int i=0 ; i<100 ; i++ )
	{
		dpm.Updata();
	}

	// 学習結果保存
	dpm.SaveModel( "result" );

	Free( data );

	return 0;
}

