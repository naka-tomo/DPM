
/*******************************************************************

DPMTest.cpp		developed by naka_t	2011.02.23

	DPMクラスの使用例

  Copyright (C) 2011  naka_t <naka_t@apple.ee.uec.ac.jp>
 *******************************************************************/
#include "stdafx.h"
#include "GaussianDPM.h"
#include "utility.h"


int _tmain(int argc, _TCHAR* argv[])
{
	CGaussianDPM dpm;			// DPMクラス
	int num, dim;				// データ数と次元
	std::vector<int> cluster;	//

	// データの読み込み
	double **data = LoadMatrix<double>( dim , num , "sample.txt" );

	// データを渡す
	dpm.SetData( data , num , dim );

	// パラメタの更新
	for(int i=0 ; i<500 ; i++ )
	{
		dpm.Updata();
	}

	// 学習結果保存
	dpm.SaveModel( "result" );
	
	// 結果を表示
	cluster = dpm.GetClusteringResult();
	for(int i=0 ; i<cluster.size() ; i++ )
	{
		printf("%d	->	%d\n" , i , cluster[i] );
	}

	Free( data );

	return 0;
}

