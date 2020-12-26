//////////////////////////////////////////////////////////////////////
// Matrix.cpp
//
// ����������� CMatrix ��ʵ���ļ�
//
//////////////////////////////////////////////////////////////////////

#include "Matrix.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// �������캯��
//////////////////////////////////////////////////////////////////////
CMatrix::CMatrix()
{
	m_nNumColumns = 1;
	m_nNumRows = 1;
	m_pData = NULL;
	BOOL bSuccess = Init(m_nNumRows, m_nNumColumns);
	ASSERT(bSuccess);
}

//////////////////////////////////////////////////////////////////////
// ָ�����й��캯��
//
// ������
// 1. int nRows - ָ���ľ�������
// 2. int nCols - ָ���ľ�������
//////////////////////////////////////////////////////////////////////
CMatrix::CMatrix(int nRows, int nCols)
{
	m_nNumRows = nRows;
	m_nNumColumns = nCols;
	m_pData = NULL;
	BOOL bSuccess = Init(m_nNumRows, m_nNumColumns);
	ASSERT(bSuccess);
}

//////////////////////////////////////////////////////////////////////
// ��ʼ������
//
// ������
// 1. int nRows - ָ���ľ�������
// 2. int nCols - ָ���ľ�������
//
// ����ֵ��BOOL �ͣ���ʼ���Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////

BOOL CMatrix::Init(int nRows, int nCols)
{
	if (m_pData)
	{
		delete[] m_pData;
		m_pData = NULL;
	}

	m_nNumRows = nRows;
	m_nNumColumns = nCols;
	int nSize = nCols*nRows;
	if (nSize < 0)
		return FALSE;

	// �����ڴ�
	m_pData = new double[nSize];
	
	if (m_pData == NULL)
		return FALSE;					// �ڴ����ʧ��
	if (IsBadReadPtr(m_pData, sizeof(double) * nSize))
		return FALSE;

	// ����Ԫ��ֵ��0
		
	memset(m_pData, 0, sizeof(double) * nSize);

	return TRUE;
}
//////////////////////////////////////////////////////////////////////
// ָ��ֵ���캯��
//
// ������
// 1. int nRows - ָ���ľ�������
// 2. int nCols - ָ���ľ�������
// 3. double value[] - һά���飬����ΪnRows*nCols���洢�����Ԫ�ص�ֵ
//////////////////////////////////////////////////////////////////////
CMatrix::CMatrix(int nRows, int nCols, double value[])
{
	m_nNumRows = nRows;
	m_nNumColumns = nCols;
	m_pData = NULL;
	BOOL bSuccess = Init(m_nNumRows, m_nNumColumns);
	ASSERT(bSuccess);

    /********************************************************/

	SetData(value);

	/********************************************************/

}

//////////////////////////////////////////////////////////////////////
// ���þ����Ԫ�ص�ֵ
//
// ������
// 1. double value[] - һά���飬����Ϊm_nNumColumns*m_nNumRows���洢
//                     �����Ԫ�ص�ֵ
//
// ����ֵ����
//////////////////////////////////////////////////////////////////////
void CMatrix::SetData(double value[])
{
	// empty the memory
	memset(m_pData, 0, sizeof(double) * m_nNumColumns*m_nNumRows);
	// copy data
	memcpy(m_pData, value, sizeof(double)*m_nNumColumns*m_nNumRows);
}

//////////////////////////////////////////////////////////////////////
// �����캯��
//
// ������
// 1. int nSize - ����������
//////////////////////////////////////////////////////////////////////
CMatrix::CMatrix(int nSize)
{
	m_nNumRows = nSize;
	m_nNumColumns = nSize;
	m_pData = NULL;
	BOOL bSuccess = Init(nSize, nSize);
	ASSERT (bSuccess);
}

//////////////////////////////////////////////////////////////////////
// �����캯��
//
// ������
// 1. int nSize - ����������
// 2. double value[] - һά���飬����ΪnRows*nRows���洢�����Ԫ�ص�ֵ
//////////////////////////////////////////////////////////////////////
CMatrix::CMatrix(int nSize, double value[])
{
	m_nNumRows = nSize;
	m_nNumColumns = nSize;
	m_pData = NULL;
	BOOL bSuccess = Init(nSize, nSize);
	ASSERT (bSuccess);

	SetData(value);
}

//////////////////////////////////////////////////////////////////////
// �������캯��
//
// ������
// 1. const CMatrix& other - Դ����
//////////////////////////////////////////////////////////////////////
CMatrix::CMatrix(const CMatrix& other)
{
	m_nNumColumns = other.GetNumColumns();
	m_nNumRows = other.GetNumRows();
	m_pData = NULL;
	BOOL bSuccess = Init(m_nNumRows, m_nNumColumns);
	ASSERT(bSuccess);

	// copy the pointer
	memcpy(m_pData, other.m_pData, sizeof(double)*m_nNumColumns*m_nNumRows);
}

//////////////////////////////////////////////////////////////////////
// ��������
//////////////////////////////////////////////////////////////////////
CMatrix::~CMatrix()
{
	if (m_pData)
	{
		delete[] m_pData;
		m_pData = NULL;
	}
}



//////////////////////////////////////////////////////////////////////
// �������ʼ��Ϊ��λ����
//
// ������
// 1. int nSize - ����������
//
// ����ֵ��BOOL �ͣ���ʼ���Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::MakeUnitMatrix(int nSize)
{
	if (! Init(nSize, nSize))
		return FALSE;

	for (int i=0; i<nSize; ++i)
		for (int j=0; j<nSize; ++j)
			if (i == j)
				SetElement(i, j, 1);

	return TRUE;
}

//////////////////////////////////////////////////////////////////////
// ���ַ���ת��Ϊ�����ֵ
//
// ������
// 1. CString s - ���ֺͷָ������ɵ��ַ���
// 2. const CString& sDelim - ����֮��ķָ�����Ĭ��Ϊ�ո�
// 3. BOOL bLineBreak - ������֮���Ƿ��лس����з���Ĭ��Ϊ��(�л��з�)
//         ���ò���ΪFALSEʱ������Ԫ��ֵ����һ�������룬�ַ����ĵ�һ��
//         ��ֵӦΪ������������ڶ�����ֵӦΪ���������
//
// ����ֵ��BOOL �ͣ�ת���Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::FromString(CString s, const CString& sDelim /*= " "*/, BOOL bLineBreak /*= TRUE*/)
{
	if (s.IsEmpty())
		return FALSE;

	// ���д���
	if (bLineBreak)
	{
		CTokenizer tk(s, "\r\n");

		CStringList ListRow;
		CString sRow;
		while (tk.Next(sRow))
		{
			sRow.TrimLeft();     //��ȥ��ߵĿո�
			sRow.TrimRight();    //��ȥ�ұߵĿո�
			if (sRow.IsEmpty())
				break;

			ListRow.AddTail(sRow);
		}

		// ����
		m_nNumRows = ListRow.GetCount();

		sRow = ListRow.GetHead();
		CTokenizer tkRow(sRow, sDelim);
		CString sElement;
		// ����
		m_nNumColumns = 0;
		while (tkRow.Next(sElement))
		{
			m_nNumColumns++;
		}

		// ��ʼ������
		if (! Init(m_nNumRows, m_nNumColumns))
			return FALSE;

		// ����ֵ
		POSITION pos = ListRow.GetHeadPosition();
		for (int i=0; i<m_nNumRows; i++)
		{
			sRow = ListRow.GetNext(pos);
			int j = 0;
			CTokenizer tkRow(sRow, sDelim);
			while (tkRow.Next(sElement))
			{
				sElement.TrimLeft();
				sElement.TrimRight();
				double v = atof(sElement);
				SetElement(i, j++, v);
			}
		}

		return TRUE;
	}
	
	// �����У����У�����

	CTokenizer tk(s, sDelim);

	CString sElement;
	
	// ����
	tk.Next(sElement);
	sElement.TrimLeft();
	sElement.TrimRight();
	m_nNumRows = atoi(sElement);

	// ����
	tk.Next(sElement);
	sElement.TrimLeft();
	sElement.TrimRight();
	m_nNumColumns = atoi(sElement);

	// ��ʼ������
	if (! Init(m_nNumRows, m_nNumColumns))
		return FALSE;

	// ����ֵ
	int i = 0, j = 0;
	while (tk.Next(sElement))
	{
		sElement.TrimLeft();
		sElement.TrimRight();
		double v = atof(sElement);
		SetElement(i, j++, v);

		if (j == m_nNumColumns)
		{
			j = 0;
			i++;
			if (i == m_nNumRows)
				break;
		}
	}

	return TRUE;
}

//////////////////////////////////////////////////////////////////////
// �������Ԫ�ص�ֵת��Ϊ�ַ���
//
// ������
// 1. const CString& sDelim - ����֮��ķָ�����Ĭ��Ϊ�ո�
// 2 BOOL bLineBreak - ������֮���Ƿ��лس����з���Ĭ��Ϊ��(�л��з�)
//
// ����ֵ��CString �ͣ�ת���õ����ַ���
//////////////////////////////////////////////////////////////////////
CString CMatrix::ToString(const CString& sDelim /*= " "*/, BOOL bLineBreak /*= TRUE*/) const
{
	CString s="";

	for (int i=0; i<m_nNumRows; ++i)
	{
		for (int j=0; j<m_nNumColumns; ++j)
		{
			CString ss;
			ss.Format("%0.8lf", GetElement(i, j));
			s += ss;

			if (bLineBreak)
			{
				if (j != m_nNumColumns-1)
					s += sDelim;
			}
			else
			{
				if (i != m_nNumRows-1 || j != m_nNumColumns-1)
					s += sDelim;
			}
		}
		if (bLineBreak)
			if (i != m_nNumRows-1)
				s += "\r\n";
	}

	return s;
}

//////////////////////////////////////////////////////////////////////
// ������ָ�����и�Ԫ�ص�ֵת��Ϊ�ַ���
//
// ������
// 1. int nRow - ָ���ľ����У�nRow = 0��ʾ��һ��
// 2. const CString& sDelim - ����֮��ķָ�����Ĭ��Ϊ�ո�
//
// ����ֵ��CString �ͣ�ת���õ����ַ���
//////////////////////////////////////////////////////////////////////
CString CMatrix::RowToString(int nRow, const CString& sDelim /*= " "*/) const
{
	CString s = "";

	if (nRow >= m_nNumRows)
		return s;

	for (int j=0; j<m_nNumColumns; ++j)
	{
		CString ss;
		ss.Format("%lf", GetElement(nRow, j));
		s += ss;
		if (j != m_nNumColumns-1)
			s += sDelim;
	}

	return s;
}

//////////////////////////////////////////////////////////////////////
// ������ָ�����и�Ԫ�ص�ֵת��Ϊ�ַ���
//
// ������
// 1. int nCol - ָ���ľ����У�nCol = 0��ʾ��һ��
// 2. const CString& sDelim - ����֮��ķָ�����Ĭ��Ϊ�ո�
//
// ����ֵ��CString �ͣ�ת���õ����ַ���
//////////////////////////////////////////////////////////////////////
CString CMatrix::ColToString(int nCol, const CString& sDelim /*= " "*/) const
{
	CString s = "";

	if (nCol >= m_nNumColumns)
		return s;

	for (int i=0; i<m_nNumRows; ++i)
	{
		CString ss;
		ss.Format("%lf", GetElement(i, nCol));
		s += ss;
		if (i != m_nNumRows-1)
			s += sDelim;
	}

	return s;
}



//////////////////////////////////////////////////////////////////////
// ����ָ��Ԫ�ص�ֵ
//
// ������
// 1. int nRows - ָ���ľ�������
// 2. int nCols - ָ���ľ�������
// 3. double value - ָ��Ԫ�ص�ֵ
//
// ����ֵ��BOOL �ͣ�˵�������Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::SetElement(int nRow, int nCol, double value)
{
	if (nCol < 0 || nCol >= m_nNumColumns || nRow < 0 || nRow >= m_nNumRows)
		return FALSE;						// array bounds error
	if (m_pData == NULL)
		return FALSE;							// bad pointer error
	
	m_pData[nCol + nRow * m_nNumColumns] = value;

	return TRUE;
}

//////////////////////////////////////////////////////////////////////
// ����ָ��Ԫ�ص�ֵ
//
// ������
// 1. int nRows - ָ���ľ�������
// 2. int nCols - ָ���ľ�������
//
// ����ֵ��double �ͣ�ָ��Ԫ�ص�ֵ
//////////////////////////////////////////////////////////////////////
double CMatrix::GetElement(int nRow, int nCol) const
{
	ASSERT(nCol >= 0 && nCol < m_nNumColumns && nRow >= 0 && nRow < m_nNumRows); // array bounds error
	ASSERT(m_pData);							// bad pointer error
	return m_pData[nCol + nRow * m_nNumColumns] ;
}

//////////////////////////////////////////////////////////////////////
// ��ȡ���������
//
// ��������
//
// ����ֵ��int �ͣ����������
//////////////////////////////////////////////////////////////////////
int	CMatrix::GetNumColumns() const
{
	return m_nNumColumns;
}

//////////////////////////////////////////////////////////////////////
// ��ȡ���������
//
// ��������
//
// ����ֵ��int �ͣ����������
//////////////////////////////////////////////////////////////////////
int	CMatrix::GetNumRows() const
{
	return m_nNumRows;
}

//////////////////////////////////////////////////////////////////////
// ��ȡ���������
//
// ��������
//
// ����ֵ��double��ָ�룬ָ������Ԫ�ص����ݻ�����
//////////////////////////////////////////////////////////////////////
double* CMatrix::GetData() const
{
	return m_pData;
}

//////////////////////////////////////////////////////////////////////
// ��ȡָ���е�����
//
// ������
// 1. int nRows - ָ���ľ�������
// 2.  double* pVector - ָ�������и�Ԫ�صĻ�����
//
// ����ֵ��int �ͣ�������Ԫ�صĸ����������������
//////////////////////////////////////////////////////////////////////
int CMatrix::GetRowVector(int nRow, double* pVector) const
{
	if (pVector == NULL)
		delete pVector;

	pVector = new double[m_nNumColumns];
	ASSERT(pVector != NULL);

	for (int j=0; j<m_nNumColumns; ++j)
		pVector[j] = GetElement(nRow, j);

	return m_nNumColumns;
}

//////////////////////////////////////////////////////////////////////
// ��ȡָ���е�����
//
// ������
// 1. int nCols - ָ���ľ�������
// 2.  double* pVector - ָ�������и�Ԫ�صĻ�����
//
// ����ֵ��int �ͣ�������Ԫ�صĸ����������������
//////////////////////////////////////////////////////////////////////
int CMatrix::GetColVector(int nCol, double* pVector) const
{
	if (pVector == NULL)
		delete pVector;

	pVector = new double[m_nNumRows];
	ASSERT(pVector != NULL);

	for (int i=0; i<m_nNumRows; ++i)
		pVector[i] = GetElement(i, nCol);

	return m_nNumRows;
}

//////////////////////////////////////////////////////////////////////
// ���������=��������ֵ
//
// ������
// 1. const CMatrix& other - ���ڸ�����ֵ��Դ����
//
// ����ֵ��CMatrix�͵����ã������õľ�����other���
//////////////////////////////////////////////////////////////////////
CMatrix& CMatrix::operator=(const CMatrix& other)
{
	if (&other != this)
	{
		BOOL bSuccess = Init(other.GetNumRows(), other.GetNumColumns());
		ASSERT(bSuccess);

		// copy the pointer
		memcpy(m_pData, other.m_pData, sizeof(double)*m_nNumColumns*m_nNumRows);
	}

	// finally return a reference to ourselves
	return *this ;
}

//////////////////////////////////////////////////////////////////////
// ���������==���жϾ����Ƿ����
//
// ������
// 1. const CMatrix& other - ���ڱȽϵľ���
//
// ����ֵ��BOOL �ͣ��������������ΪTRUE������ΪFALSE
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::operator==(const CMatrix& other) const
{
	// ���ȼ���������Ƿ����
	if (m_nNumColumns != other.GetNumColumns() || m_nNumRows != other.GetNumRows())
		return FALSE;

	for (int i=0; i<m_nNumRows; ++i)
	{
		for (int j=0; j<m_nNumColumns; ++j)
		{
			if (GetElement(i, j) != other.GetElement(i, j))
				return FALSE;
		}
	}

	return TRUE;
}

//////////////////////////////////////////////////////////////////////
// ���������!=���жϾ����Ƿ����
//
// ������
// 1. const CMatrix& other - ���ڱȽϵľ���
//
// ����ֵ��BOOL �ͣ����������������ΪTRUE������ΪFALSE
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::operator!=(const CMatrix& other) const
{
	return !(*this == other);
}

//////////////////////////////////////////////////////////////////////
// ���������+��ʵ�־���ļӷ�
//
// ������
// 1. const CMatrix& other - ��ָ��������ӵľ���
//
// ����ֵ��CMatrix�ͣ�ָ��������other���֮��
//////////////////////////////////////////////////////////////////////
CMatrix	CMatrix::operator+(const CMatrix& other) const
{
	// ���ȼ���������Ƿ����
	ASSERT (m_nNumColumns == other.GetNumColumns() && m_nNumRows == other.GetNumRows());

	// ����������
	CMatrix	result(*this) ;		// ��������
	// ����ӷ�
	for (int i = 0 ; i < m_nNumRows ; ++i)
	{
		for (int j = 0 ; j <  m_nNumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) + other.GetElement(i, j)) ;
	}

	return result ;
}

//////////////////////////////////////////////////////////////////////
// ���������-��ʵ�־���ļ���
//
// ������
// 1. const CMatrix& other - ��ָ����������ľ���
//
// ����ֵ��CMatrix�ͣ�ָ��������other���֮��
//////////////////////////////////////////////////////////////////////
CMatrix	CMatrix::operator-(const CMatrix& other) const
{
	// ���ȼ���������Ƿ����
	ASSERT (m_nNumColumns == other.GetNumColumns() && m_nNumRows == other.GetNumRows());

	// ����Ŀ�����
	CMatrix	result(*this) ;		// copy ourselves
	// ���м�������
	for (int i = 0 ; i < m_nNumRows ; ++i)
	{
		for (int j = 0 ; j <  m_nNumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) - other.GetElement(i, j)) ;
	}

	return result ;
}

//////////////////////////////////////////////////////////////////////
// ���������*��ʵ�־��������
//
// ������
// 1. double value - ��ָ��������˵�ʵ��
//
// ����ֵ��CMatrix�ͣ�ָ��������value���֮��
//////////////////////////////////////////////////////////////////////
CMatrix	CMatrix::operator*(double value) const
{
	// ����Ŀ�����
	CMatrix	result(*this) ;		// copy ourselves
	// ��������
	for (int i = 0 ; i < m_nNumRows ; ++i)
	{
		for (int j = 0 ; j <  m_nNumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) * value) ;
	}

	return result ;
}

//////////////////////////////////////////////////////////////////////
// ���������*��ʵ�־���ĳ˷�
//
// ������
// 1. const CMatrix& other - ��ָ��������˵ľ���
//
// ����ֵ��CMatrix�ͣ�ָ��������other���֮��
//////////////////////////////////////////////////////////////////////
CMatrix	CMatrix::operator*(const CMatrix& other) const
{
	// ���ȼ���������Ƿ����Ҫ��
	ASSERT (m_nNumColumns == other.GetNumRows());

	// construct the object we are going to return
	CMatrix	result(m_nNumRows, other.GetNumColumns()) ;

	// ����˷�����
	//
	// [A][B][C]   [G][H]     [A*G + B*I + C*K][A*H + B*J + C*L]
	// [D][E][F] * [I][J] =   [D*G + E*I + F*K][D*H + E*J + F*L]
	//             [K][L]
	//
	double	value ;
	for (int i = 0 ; i < result.GetNumRows() ; ++i)
	{
		for (int j = 0 ; j < other.GetNumColumns() ; ++j)
		{
			value = 0.0 ;
			for (int k = 0 ; k < m_nNumColumns ; ++k)
			{
				value += GetElement(i, k) * other.GetElement(k, j) ;
			}

			result.SetElement(i, j, value) ;
		}
	}

	return result ;
}

//////////////////////////////////////////////////////////////////////
// �����ת��
//
// ��������
//
// ����ֵ��CMatrix�ͣ�ָ������ת�þ���
//////////////////////////////////////////////////////////////////////
CMatrix CMatrix::Transpose() const
{
	// ����Ŀ�����
	CMatrix	Trans(m_nNumColumns, m_nNumRows);

	// ת�ø�Ԫ��
	for (int i = 0 ; i < m_nNumRows ; ++i)
	{
		for (int j = 0 ; j < m_nNumColumns ; ++j)
			Trans.SetElement(j, i, GetElement(i, j)) ;
	}

	return Trans;
}

//////////////////////////////////////////////////////////////////////
// ʵ���������ȫѡ��Ԫ��˹��Լ����
//
// ��������
//
// ����ֵ��BOOL�ͣ������Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::InvertGaussJordan()
{
	int *pnRow, *pnCol,i,j,k,l,u,v;
    double d = 0, p = 0;

	// �����ڴ�
    pnRow = new int[m_nNumColumns];
    pnCol = new int[m_nNumColumns];
	if (pnRow == NULL || pnCol == NULL)
		return FALSE;

	// ��Ԫ
    for (k=0; k<=m_nNumColumns-1; k++)
    { 
		d=0.0;
        for (i=k; i<=m_nNumColumns-1; i++)
		{
			for (j=k; j<=m_nNumColumns-1; j++)
			{ 
				l=i*m_nNumColumns+j; p=fabs(m_pData[l]);
				if (p>d) 
				{ 
					d=p; 
					pnRow[k]=i; 
					pnCol[k]=j;
				}
			}
		}
        
		// ʧ��
		if (d == 0.0)
		{
			delete[] pnRow;
			delete[] pnCol;
			return FALSE;
		}

        if (pnRow[k] != k)
		{
			for (j=0; j<=m_nNumColumns-1; j++)
			{ 
				u=k*m_nNumColumns+j; 
				v=pnRow[k]*m_nNumColumns+j;
				p=m_pData[u]; 
				m_pData[u]=m_pData[v]; 
				m_pData[v]=p;
			}
		}
        
		if (pnCol[k] != k)
		{
			for (i=0; i<=m_nNumColumns-1; i++)
            { 
				u=i*m_nNumColumns+k; 
				v=i*m_nNumColumns+pnCol[k];
				p=m_pData[u]; 
				m_pData[u]=m_pData[v]; 
				m_pData[v]=p;
            }
		}

        l=k*m_nNumColumns+k;
        m_pData[l]=1.0/m_pData[l];
        for (j=0; j<=m_nNumColumns-1; j++)
		{
			if (j != k)
            { 
				u=k*m_nNumColumns+j; 
				m_pData[u]=m_pData[u]*m_pData[l];
			}
		}

        for (i=0; i<=m_nNumColumns-1; i++)
		{
			if (i!=k)
			{
				for (j=0; j<=m_nNumColumns-1; j++)
				{
					if (j!=k)
					{ 
						u=i*m_nNumColumns+j;
						m_pData[u]=m_pData[u]-m_pData[i*m_nNumColumns+k]*m_pData[k*m_nNumColumns+j];
					}
                }
			}
		}

        for (i=0; i<=m_nNumColumns-1; i++)
		{
			if (i!=k)
            { 
				u=i*m_nNumColumns+k; 
				m_pData[u]=-m_pData[u]*m_pData[l];
			}
		}
    }

    // �����ָ����д���
    for (k=m_nNumColumns-1; k>=0; k--)
    { 
		if (pnCol[k]!=k)
		{
			for (j=0; j<=m_nNumColumns-1; j++)
            { 
				u=k*m_nNumColumns+j; 
				v=pnCol[k]*m_nNumColumns+j;
				p=m_pData[u]; 
				m_pData[u]=m_pData[v]; 
				m_pData[v]=p;
            }
		}

        if (pnRow[k]!=k)
		{
			for (i=0; i<=m_nNumColumns-1; i++)
            { 
				u=i*m_nNumColumns+k; 
				v=i*m_nNumColumns+pnRow[k];
				p=m_pData[u]; 
				m_pData[u]=m_pData[v]; 
				m_pData[v]=p;
            }
		}
    }

	// �����ڴ�
	delete[] pnRow;
	delete[] pnCol;

	// �ɹ�����
	return TRUE;
}


//////////////////////////////////////////////////////////////////////
// �Գ��������������
//
// ��������
//
// ����ֵ��BOOL�ͣ������Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::InvertSsgj()
{ 
	int i, j ,k, m;
    double w, g, *pTmp;

	// ��ʱ�ڴ�
    pTmp = new double[m_nNumColumns];

	// ���д���
    for (k=0; k<=m_nNumColumns-1; k++)
    { 
		w=m_pData[0];
        if (w == 0.0)
        { 
			delete[] pTmp;
			return FALSE;
		}

        m=m_nNumColumns-k-1;
        for (i=1; i<=m_nNumColumns-1; i++)
        { 
			g=m_pData[i*m_nNumColumns]; 
			pTmp[i]=g/w;
            if (i<=m) 
				pTmp[i]=-pTmp[i];
            for (j=1; j<=i; j++)
              m_pData[(i-1)*m_nNumColumns+j-1]=m_pData[i*m_nNumColumns+j]+g*pTmp[j];
        }

        m_pData[m_nNumColumns*m_nNumColumns-1]=1.0/w;
        for (i=1; i<=m_nNumColumns-1; i++)
			m_pData[(m_nNumColumns-1)*m_nNumColumns+i-1]=pTmp[i];
    }

	// ���е���
    for (i=0; i<=m_nNumColumns-2; i++)
		for (j=i+1; j<=m_nNumColumns-1; j++)
			m_pData[i*m_nNumColumns+j]=m_pData[j*m_nNumColumns+i];

	// ��ʱ�ڴ�����
	delete[] pTmp;

	return TRUE;
}

              










