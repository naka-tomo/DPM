#include "StdAfx.h"
#include "GaussianTable.h"
#include "gamma.h"
#include <math.h>

GaussianParams::GaussianParams(int dim_)
{
	Create( dim_ );
}

void GaussianParams::Init( double **data , int num )
{
	r = 1;
	nu = dim + 2;
	N = 0;

	// データ全体の平均を初期値に
	for(int n=0 ; n<num ; n++ )
	{
		cv::Mat x( dim , 1 , CV_64F , data[n] );
		m += x;
	}
	m /= num;


	// データの分散を初期値に
	for(int n=0 ; n<num ; n++ )
	{
		cv::Mat x( dim , 1 , CV_64F , data[n] );
		S += (x-m) * ((x-m).t());		
	}
	S /= num;


	// |S| = 0 となると計算できないので初期値を変更
	if( cv::determinant(S) <= 0.000001 )
	{
		// 対角成分だけを残す(これで良いかは要検討・・・)
		for(int d1=0 ; d1<dim ; d1++ )
			for(int d2=0 ; d2<dim ; d2++ )
				if( d1 != d2 ) S.at<double>(d1,d2) = 0;

		if( cv::determinant(S) <= 0 )
		{
			MessageBox( 0 , "初期値の設定に失敗。計算出来ません。データのスケールを変更することで計算できる可能性があります。" , "Error" , 0 );
		}
	}
}

void GaussianParams::Create( int dim_ )
{
	dim = dim_;
	X = cv::Mat( dim , 1 , CV_64F , cv::Scalar(0) );
	C = cv::Mat( dim , dim , CV_64F , cv::Scalar(0) );
	m = cv::Mat( dim , 1 , CV_64F , cv::Scalar(0) );
	S = cv::Mat( dim , dim , CV_64F , cv::Scalar(0) );
	r = 0;
	nu = 0;
	N = 0;
}


GaussianParams GaussianParams::Clone()
{
	GaussianParams p;
	p.X = X.clone();
	p.C = C.clone();
	p.m = m.clone();
	p.S = S.clone();
	p.r = r;
	p.nu = nu;
	p.N = N;
	p.dim = dim;
	return p;
}

CGaussianTable::CGaussianTable( GaussianParams &init )
{
	m_init = init.Clone();	
	m_param = init.Clone();
}

CGaussianTable::~CGaussianTable(void)
{
}


void CGaussianTable::AddData( double *data )
{
	AddData( m_param , data );
}

void CGaussianTable::AddData( GaussianParams &p , double *data )
{
	cv::Mat x( m_param.dim , 1 , CV_64F , data );

	p.X += x;
	p.C += x*x.t();
	p.r++;
	p.nu++;
	p.N++;
	p.m = (p.X + m_init.r * m_init.m)/(m_init.r + p.N );
	p.S = - p.r * p.m * p.m.t() +
				p.C + 
				m_init.S + 
				m_init.r * m_init.m * m_init.m.t();

}

void CGaussianTable::DeleteData( double *data )
{
	cv::Mat x( m_param.dim , 1 , CV_64F , data );

	m_param.X -= x;
	m_param.C -= x*x.t();
	m_param.r--;
	m_param.nu--;
	m_param.N--;
	m_param.m = (m_param.X + m_init.r * m_init.m)/(m_init.r + m_param.N );
	m_param.S = - m_param.r * m_param.m * m_param.m.t() +
				m_param.C + m_init.S + 
				m_init.r * m_init.m * m_init.m.t();
}


double CGaussianTable::CalcLogPx( double *addData )
{
	GaussianParams p;
	double P = 0;

	if( addData ) 
	{
		p = m_param.Clone();	// 複製
		AddData( p , addData );	// データを追加
	}
	else
	{
		p = m_param;			// パラメタを参照
	}

	P = - p.N * p.dim * 0.5 * log( M_PI ) - 
		p.dim * 0.5 * log( (double)p.r) -
		p.nu * 0.5 * log( cv::determinant( p.S ) );

	if( !_finite( P ) )
	{
		printf("%lf\n" , cv::determinant( p.S ) );
		MessageBox( 0 , "計算出来ません。データのスケールを変更することで計算できる可能性があります。" , "Error" , 0 );
	}


	for(int d=1 ; d<=p.dim ; d++ )
	{
		P += loggamma( 0.5*(p.nu+1-d) );
	}

	return P;
}

double CGaussianTable::CalcLogLikilihood( double *data )
{
	// log( P(x_new,X) ) - log( P(X) )
	return CalcLogPx( data ) - CalcLogPx( NULL );
}