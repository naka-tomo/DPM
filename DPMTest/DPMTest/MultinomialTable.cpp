#include "StdAfx.h"
#include "MultinomialTable.h"
#include "gamma.h"
#include <math.h>

CMultinomialTable::CMultinomialTable(int dim) 
: m_alpha( dim , DPM_ALPHA ) , m_numData(0) , m_dim(dim)
{
}

CMultinomialTable::~CMultinomialTable(void)
{
}


void CMultinomialTable::AddData( double *data )
{
	m_numData ++;

	double *p = &m_alpha[0];
	for(int i=0 ; i<m_dim ; i++ ) p[i] += data[i];
}

void CMultinomialTable::DeleteData( double *data )
{
	m_numData --;

	double *p = &m_alpha[0];
	for(int i=0 ; i<m_dim ; i++ ) p[i] -= data[i];
}

double CMultinomialTable::CalcLogZ( double *alpha , double *addData )
{
	double sum = 0;
	double z = 0;

	for(int i=0 ; i<m_dim ; i++ )
	{
		if( addData )
		{
			sum += (alpha[i] + addData[i]); 
			z += loggamma( alpha[i] + addData[i] );
		}
		else
		{
			sum += alpha[i];
			z += loggamma( alpha[i] );
		}

	}

	z -= loggamma( sum );

	return z;
}

double CMultinomialTable::CalcLogLikilihood( double *data )
{
	return CalcLogZ( &m_alpha[0] , data ) - CalcLogZ( &m_alpha[0] );
}