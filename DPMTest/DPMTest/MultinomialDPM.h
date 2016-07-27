#pragma once
/*******************************************************************

MultinomialDPM.h		developed by naka_t	2011.02.01

	多項分布 Dirichlet Process Mixture

　＊公開用にプログラムを整理				naka_t	2010.02.01

  Copyright (C) 2011  naka_t <naka_t@apple.ee.uec.ac.jp>
 *******************************************************************/
#include "MultinomialTable.h"

#define MAX_TABLES 100				// テーブルの最大数（計算用バッファ確保に必要）
#define DP_CONCPARA_PRIOR_A 0.1		// Concentration parameterのガンマ事前分布のパラメタ
#define DP_CONCPARA_PRIOR_B 0.1		// Concentration parameterのガンマ事前分布のパラメタ

class CMultinomialDPM
{
public:
	CMultinomialDPM(void);
	~CMultinomialDPM(void);

	void SetData( double **data , int numData , int dim );	// データのセット
	void Updata();											// パラメータ更新
	void SaveModel( const char *dir );						// モデルの保存

protected:
	void Release();								// メモリ解放
	int Sampling( double *data );				// 客が座るテーブルのサンプリング
	double SamplingGamma( double oldGamma );	// concentrate paramter のサンプリング
	void DeleteEmptyTables();					// 空のテーブルを削除
	void DeleteTable( int t );					// 特定のテーブルを削除
	double Rand();								// 0～1の乱数発生

	int m_numData;								// データ数
	int m_dimData;								// データの次元
	double **m_data;							// データ
	int *m_tableID;								// データが座ったテーブル
	double m_gamma;
	std::vector<CMultinomialTable> m_tables;	// データが座るテーブル
	std::vector<int> m_numTables;				// 各回でのテーブルの数
	double m_Pz[MAX_TABLES];					// 計算用バッファ

};

