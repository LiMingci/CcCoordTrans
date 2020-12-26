//////////////////////////////////////////////////////////////////////
// Matrix.cpp
//
// 操作矩阵的类 CMatrix 的实现文件
//
//////////////////////////////////////////////////////////////////////

#include "Matrix.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// 基本构造函数
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
// 指定行列构造函数
//
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
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
// 初始化函数
//
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
//
// 返回值：BOOL 型，初始化是否成功
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

	// 分配内存
	m_pData = new double[nSize];
	
	if (m_pData == NULL)
		return FALSE;					// 内存分配失败
	if (IsBadReadPtr(m_pData, sizeof(double) * nSize))
		return FALSE;

	// 将各元素值置0
		
	memset(m_pData, 0, sizeof(double) * nSize);

	return TRUE;
}
//////////////////////////////////////////////////////////////////////
// 指定值构造函数
//
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
// 3. double value[] - 一维数组，长度为nRows*nCols，存储矩阵各元素的值
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
// 设置矩阵各元素的值
//
// 参数：
// 1. double value[] - 一维数组，长度为m_nNumColumns*m_nNumRows，存储
//                     矩阵各元素的值
//
// 返回值：无
//////////////////////////////////////////////////////////////////////
void CMatrix::SetData(double value[])
{
	// empty the memory
	memset(m_pData, 0, sizeof(double) * m_nNumColumns*m_nNumRows);
	// copy data
	memcpy(m_pData, value, sizeof(double)*m_nNumColumns*m_nNumRows);
}

//////////////////////////////////////////////////////////////////////
// 方阵构造函数
//
// 参数：
// 1. int nSize - 方阵行列数
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
// 方阵构造函数
//
// 参数：
// 1. int nSize - 方阵行列数
// 2. double value[] - 一维数组，长度为nRows*nRows，存储方阵各元素的值
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
// 拷贝构造函数
//
// 参数：
// 1. const CMatrix& other - 源矩阵
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
// 析构函数
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
// 将方阵初始化为单位矩阵
//
// 参数：
// 1. int nSize - 方阵行列数
//
// 返回值：BOOL 型，初始化是否成功
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
// 将字符串转化为矩阵的值
//
// 参数：
// 1. CString s - 数字和分隔符构成的字符串
// 2. const CString& sDelim - 数字之间的分隔符，默认为空格
// 3. BOOL bLineBreak - 行与行之间是否有回车换行符，默认为真(有换行符)
//         当该参数为FALSE时，所有元素值都在一行中输入，字符串的第一个
//         数值应为矩阵的行数，第二个数值应为矩阵的列数
//
// 返回值：BOOL 型，转换是否成功
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::FromString(CString s, const CString& sDelim /*= " "*/, BOOL bLineBreak /*= TRUE*/)
{
	if (s.IsEmpty())
		return FALSE;

	// 分行处理
	if (bLineBreak)
	{
		CTokenizer tk(s, "\r\n");

		CStringList ListRow;
		CString sRow;
		while (tk.Next(sRow))
		{
			sRow.TrimLeft();     //除去左边的空格
			sRow.TrimRight();    //除去右边的空格
			if (sRow.IsEmpty())
				break;

			ListRow.AddTail(sRow);
		}

		// 行数
		m_nNumRows = ListRow.GetCount();

		sRow = ListRow.GetHead();
		CTokenizer tkRow(sRow, sDelim);
		CString sElement;
		// 列数
		m_nNumColumns = 0;
		while (tkRow.Next(sElement))
		{
			m_nNumColumns++;
		}

		// 初始化矩阵
		if (! Init(m_nNumRows, m_nNumColumns))
			return FALSE;

		// 设置值
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
	
	// 不分行（单行）处理

	CTokenizer tk(s, sDelim);

	CString sElement;
	
	// 行数
	tk.Next(sElement);
	sElement.TrimLeft();
	sElement.TrimRight();
	m_nNumRows = atoi(sElement);

	// 列数
	tk.Next(sElement);
	sElement.TrimLeft();
	sElement.TrimRight();
	m_nNumColumns = atoi(sElement);

	// 初始化矩阵
	if (! Init(m_nNumRows, m_nNumColumns))
		return FALSE;

	// 设置值
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
// 将矩阵各元素的值转化为字符串
//
// 参数：
// 1. const CString& sDelim - 数字之间的分隔符，默认为空格
// 2 BOOL bLineBreak - 行与行之间是否有回车换行符，默认为真(有换行符)
//
// 返回值：CString 型，转换得到的字符串
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
// 将矩阵指定行中各元素的值转化为字符串
//
// 参数：
// 1. int nRow - 指定的矩阵行，nRow = 0表示第一行
// 2. const CString& sDelim - 数字之间的分隔符，默认为空格
//
// 返回值：CString 型，转换得到的字符串
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
// 将矩阵指定列中各元素的值转化为字符串
//
// 参数：
// 1. int nCol - 指定的矩阵行，nCol = 0表示第一列
// 2. const CString& sDelim - 数字之间的分隔符，默认为空格
//
// 返回值：CString 型，转换得到的字符串
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
// 设置指定元素的值
//
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
// 3. double value - 指定元素的值
//
// 返回值：BOOL 型，说明设置是否成功
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
// 设置指定元素的值
//
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
//
// 返回值：double 型，指定元素的值
//////////////////////////////////////////////////////////////////////
double CMatrix::GetElement(int nRow, int nCol) const
{
	ASSERT(nCol >= 0 && nCol < m_nNumColumns && nRow >= 0 && nRow < m_nNumRows); // array bounds error
	ASSERT(m_pData);							// bad pointer error
	return m_pData[nCol + nRow * m_nNumColumns] ;
}

//////////////////////////////////////////////////////////////////////
// 获取矩阵的列数
//
// 参数：无
//
// 返回值：int 型，矩阵的列数
//////////////////////////////////////////////////////////////////////
int	CMatrix::GetNumColumns() const
{
	return m_nNumColumns;
}

//////////////////////////////////////////////////////////////////////
// 获取矩阵的行数
//
// 参数：无
//
// 返回值：int 型，矩阵的行数
//////////////////////////////////////////////////////////////////////
int	CMatrix::GetNumRows() const
{
	return m_nNumRows;
}

//////////////////////////////////////////////////////////////////////
// 获取矩阵的数据
//
// 参数：无
//
// 返回值：double型指针，指向矩阵各元素的数据缓冲区
//////////////////////////////////////////////////////////////////////
double* CMatrix::GetData() const
{
	return m_pData;
}

//////////////////////////////////////////////////////////////////////
// 获取指定行的向量
//
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2.  double* pVector - 指向向量中各元素的缓冲区
//
// 返回值：int 型，向量中元素的个数，即矩阵的列数
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
// 获取指定列的向量
//
// 参数：
// 1. int nCols - 指定的矩阵列数
// 2.  double* pVector - 指向向量中各元素的缓冲区
//
// 返回值：int 型，向量中元素的个数，即矩阵的行数
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
// 重载运算符=，给矩阵赋值
//
// 参数：
// 1. const CMatrix& other - 用于给矩阵赋值的源矩阵
//
// 返回值：CMatrix型的引用，所引用的矩阵与other相等
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
// 重载运算符==，判断矩阵是否相等
//
// 参数：
// 1. const CMatrix& other - 用于比较的矩阵
//
// 返回值：BOOL 型，两个矩阵相等则为TRUE，否则为FALSE
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::operator==(const CMatrix& other) const
{
	// 首先检查行列数是否相等
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
// 重载运算符!=，判断矩阵是否不相等
//
// 参数：
// 1. const CMatrix& other - 用于比较的矩阵
//
// 返回值：BOOL 型，两个不矩阵相等则为TRUE，否则为FALSE
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::operator!=(const CMatrix& other) const
{
	return !(*this == other);
}

//////////////////////////////////////////////////////////////////////
// 重载运算符+，实现矩阵的加法
//
// 参数：
// 1. const CMatrix& other - 与指定矩阵相加的矩阵
//
// 返回值：CMatrix型，指定矩阵与other相加之和
//////////////////////////////////////////////////////////////////////
CMatrix	CMatrix::operator+(const CMatrix& other) const
{
	// 首先检查行列数是否相等
	ASSERT (m_nNumColumns == other.GetNumColumns() && m_nNumRows == other.GetNumRows());

	// 构造结果矩阵
	CMatrix	result(*this) ;		// 拷贝构造
	// 矩阵加法
	for (int i = 0 ; i < m_nNumRows ; ++i)
	{
		for (int j = 0 ; j <  m_nNumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) + other.GetElement(i, j)) ;
	}

	return result ;
}

//////////////////////////////////////////////////////////////////////
// 重载运算符-，实现矩阵的减法
//
// 参数：
// 1. const CMatrix& other - 与指定矩阵相减的矩阵
//
// 返回值：CMatrix型，指定矩阵与other相减之差
//////////////////////////////////////////////////////////////////////
CMatrix	CMatrix::operator-(const CMatrix& other) const
{
	// 首先检查行列数是否相等
	ASSERT (m_nNumColumns == other.GetNumColumns() && m_nNumRows == other.GetNumRows());

	// 构造目标矩阵
	CMatrix	result(*this) ;		// copy ourselves
	// 进行减法操作
	for (int i = 0 ; i < m_nNumRows ; ++i)
	{
		for (int j = 0 ; j <  m_nNumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) - other.GetElement(i, j)) ;
	}

	return result ;
}

//////////////////////////////////////////////////////////////////////
// 重载运算符*，实现矩阵的数乘
//
// 参数：
// 1. double value - 与指定矩阵相乘的实数
//
// 返回值：CMatrix型，指定矩阵与value相乘之积
//////////////////////////////////////////////////////////////////////
CMatrix	CMatrix::operator*(double value) const
{
	// 构造目标矩阵
	CMatrix	result(*this) ;		// copy ourselves
	// 进行数乘
	for (int i = 0 ; i < m_nNumRows ; ++i)
	{
		for (int j = 0 ; j <  m_nNumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) * value) ;
	}

	return result ;
}

//////////////////////////////////////////////////////////////////////
// 重载运算符*，实现矩阵的乘法
//
// 参数：
// 1. const CMatrix& other - 与指定矩阵相乘的矩阵
//
// 返回值：CMatrix型，指定矩阵与other相乘之积
//////////////////////////////////////////////////////////////////////
CMatrix	CMatrix::operator*(const CMatrix& other) const
{
	// 首先检查行列数是否符合要求
	ASSERT (m_nNumColumns == other.GetNumRows());

	// construct the object we are going to return
	CMatrix	result(m_nNumRows, other.GetNumColumns()) ;

	// 矩阵乘法，即
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
// 矩阵的转置
//
// 参数：无
//
// 返回值：CMatrix型，指定矩阵转置矩阵
//////////////////////////////////////////////////////////////////////
CMatrix CMatrix::Transpose() const
{
	// 构造目标矩阵
	CMatrix	Trans(m_nNumColumns, m_nNumRows);

	// 转置各元素
	for (int i = 0 ; i < m_nNumRows ; ++i)
	{
		for (int j = 0 ; j < m_nNumColumns ; ++j)
			Trans.SetElement(j, i, GetElement(i, j)) ;
	}

	return Trans;
}

//////////////////////////////////////////////////////////////////////
// 实矩阵求逆的全选主元高斯－约当法
//
// 参数：无
//
// 返回值：BOOL型，求逆是否成功
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::InvertGaussJordan()
{
	int *pnRow, *pnCol,i,j,k,l,u,v;
    double d = 0, p = 0;

	// 分配内存
    pnRow = new int[m_nNumColumns];
    pnCol = new int[m_nNumColumns];
	if (pnRow == NULL || pnCol == NULL)
		return FALSE;

	// 消元
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
        
		// 失败
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

    // 调整恢复行列次序
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

	// 清理内存
	delete[] pnRow;
	delete[] pnCol;

	// 成功返回
	return TRUE;
}


//////////////////////////////////////////////////////////////////////
// 对称正定矩阵的求逆
//
// 参数：无
//
// 返回值：BOOL型，求逆是否成功
//////////////////////////////////////////////////////////////////////
BOOL CMatrix::InvertSsgj()
{ 
	int i, j ,k, m;
    double w, g, *pTmp;

	// 临时内存
    pTmp = new double[m_nNumColumns];

	// 逐列处理
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

	// 行列调整
    for (i=0; i<=m_nNumColumns-2; i++)
		for (j=i+1; j<=m_nNumColumns-1; j++)
			m_pData[i*m_nNumColumns+j]=m_pData[j*m_nNumColumns+i];

	// 临时内存清理
	delete[] pTmp;

	return TRUE;
}

              










