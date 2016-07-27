#include "StdAfx.h"
#include <string.h>
#include <float.h>
#include <math.h>
#include <windows.h>
#include "randlib/randlib.h"
#include "MultinomialDPM.h"
#include "utility.h"

CMultinomialDPM::CMultinomialDPM(void)
{
	m_numData = 0;
	m_dimData = 0;
	m_data = NULL;
	m_tableID = NULL;
	m_gamma = DP_CONCPARA_PRIOR_A / DP_CONCPARA_PRIOR_B;
}

CMultinomialDPM::~CMultinomialDPM(void)
{
	Release();
}

void CMultinomialDPM::Release()
{
	delete [] m_tableID;
	m_tableID = NULL;
}

void CMultinomialDPM::SetData( double **data , int numData , int dim )
{
	Release();

	m_data = data;
	m_dimData = dim;
	m_numData = numData;
	m_tableID = new int[numData];
	for(int i=0 ; i<numData ; i++ ) m_tableID[i] = -1;

	m_tables.clear();
	m_tables.push_back( CMultinomialTable(dim) );			// 空のテーブルを追加

	m_numTables.clear();
	m_numTables.reserve(1000);
}


void CMultinomialDPM::Updata()
{
	for(int d=0 ; d<m_numData ; d++ )
	{
		int t = m_tableID[d];

		if( t>=0 ) m_tables[t].DeleteData( m_data[d] );		// データを削除
		DeleteEmptyTables();								// 空のテーブルを削除
		t = Sampling( m_data[d] );							// サンプリング
		m_tableID[d] = t;
		m_tables[t].AddData( m_data[d] );					// データ追加

		// 最後の空のテーブルに座った場合
		// 新たに空のテーブルを追加
		if( t==m_tables.size()-1 )
		{
			m_tables.push_back( CMultinomialTable(m_dimData) );
		}

		DeleteEmptyTables();					// 空のテーブルを削除
	}

	m_gamma = SamplingGamma( m_gamma );
	m_numTables.push_back( (int)m_tables.size() -1 );
}

int CMultinomialDPM::Sampling( double *data )
{
	double *Pz = m_Pz;
	int numTable = (int)m_tables.size();
	int newTable = -1;
	double max = -DBL_MAX;

	// テーブルの最後は空のテーブル
	if( numTable == 1 ) return 0;
	
	// 既存のテーブル：対数尤度
	for(int t=0 ; t<numTable ; t++ )
	{
		Pz[t] = m_tables[t].CalcLogLikilihood( data );
	}

	// 最大値を探す
	for(int t=0 ; t<numTable ; t++ ) if( max < Pz[t] ) max = Pz[t];

	// 値が小さくなりすぎるため、最大値で引く
	// 各テーブルの人気をかける
	// サンプリングのために累積確率にする
	Pz[0] = exp(Pz[0] - max) * m_tables[0].GetNumData();
	for(int t=1 ; t<numTable-1 ; t++ ) Pz[t] = Pz[t-1] + exp(Pz[t] - max) * m_tables[t].GetNumData(); 

	// 新たなテーブルを生成する確率
	Pz[numTable-1] = Pz[numTable-2] + exp(Pz[numTable-1] - max) * m_gamma;

	// サンプリングするための乱数を発生
	double rand = Rand() * Pz[numTable-1];

	// 計算した確率に従って新たなテーブルを選択
	for(newTable=0 ; newTable<numTable-1 ; newTable++ )
	{
		if( Pz[newTable] >= rand ) break;
	}

	return newTable;
}

void CMultinomialDPM::DeleteEmptyTables()
{
	int numTable = (int)m_tables.size() - 1;

	// 空のテーブルを探す
	for(int t=0 ; t<numTable ; t++ )
	{
		if( m_tables[t].GetNumData() == 0 )
		{
			DeleteTable( t );
			break;
		}
	}

}

void CMultinomialDPM::DeleteTable( int t )
{
	m_tables.erase( m_tables.begin() + t );

	// 消したテーブル分idをづらす
	for(int d=0 ; d<m_numData ; d++ )
	{
		if( m_tableID[d] > t )
		{
			m_tableID[d]--;
		}
	}
}


void CMultinomialDPM::SaveModel( const char *dir )
{
	char dirname[256];
	char filename[256];
	FILE *fpPz, *fpTable;

	strcpy( dirname , dir );
	int len = (int)strlen( dir );
	if( len==0 || dir[len-1] != '\\' || dir[len-1] != '/' ) strcat( dirname , "\\" );
	
	CreateDirectory( dirname , NULL );

	if( m_numTables.size() )
	{
		sprintf( filename , "%sNumTables.txt" , dirname );
		SaveArray( m_numTables , (int)m_numTables.size() , filename );
	}

	// 各クラスに属する確率を保存
	sprintf( filename , "%sLogLiklihood.txt" , dirname );
	fpPz = fopen( filename , "w" );
	sprintf( filename , "%sClusteringResult.txt" , dirname );
	fpTable = fopen( filename , "w" );
	for(int d=0 ; d<m_numData ; d++ )
	{
		double max = -DBL_MAX;
		int maxIdx = -1;
		for(int t=0 ; t<m_tables.size() ; t++ )
		{
			double lik = m_tables[t].CalcLogLikilihood( m_data[d] );
			fprintf( fpPz , "%lf	" ,  lik );
			if( max < lik )
			{
				max = lik;
				maxIdx = t;
			}
		}
		fprintf( fpPz , "\n" );
		fprintf( fpTable , "%d\n" , maxIdx );
	}
	fclose(fpPz);
	fclose(fpTable);
}

double CMultinomialDPM::Rand()
{
	unsigned int val;
	rand_s(&val);
	return (double)val/UINT_MAX;
}

double CMultinomialDPM::SamplingGamma( double oldGamma )
{
	for(int i=0 ; i<20 ; i++ )
	{
		float gammaB = 0;	// ガンマ関数スケールパラメータ
		float gammaA = 0;	// ガンマ関数形状パラメータ
		int numTable = (int)m_tables.size() - 1;

		// ベータ分布からサンプル生成
		float w = genbet( (float)oldGamma+1 , (float)m_numData );

		// 二値の分布をサンプリング
		int s = (Rand() * (oldGamma + numTable)) < numTable ? 1 : 0;
		gammaA = (float)(DP_CONCPARA_PRIOR_A + numTable - s);
		gammaB = (float)(DP_CONCPARA_PRIOR_B - log(w));

		// 更新
		oldGamma = (double)gengam( gammaB , gammaA );
	}

	return oldGamma;
}
