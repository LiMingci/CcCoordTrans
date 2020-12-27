//////////////////////////////////////////////////////////////////////
// Matrix.h
//
// 操作矩阵的类 CMatrix 的声明接口
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include <math.h>

#if !defined(_BITSET_)
#include <bitset>
#endif // !defined(_BITSET_)

///////////////////////////////////////////////////////////////////////////////////////
//
//(-- class CTokenizer
//
class CTokenizer
{
public:
	CTokenizer(const CString& cs, const CString& csDelim) : m_cs(cs), m_nCurPos(0)
	{
		SetDelimiters(csDelim);
	}
	void SetDelimiters(const CString& csDelim)
	{
		for (int i = 0; i < csDelim.GetLength(); ++i)
			m_sDelimiter.set(static_cast<BYTE>(csDelim[i]));
	}

	BOOL Next(CString& cs)
	{
		cs.Empty();

		while (m_nCurPos < m_cs.GetLength() && m_sDelimiter[static_cast<BYTE>(m_cs[m_nCurPos])])
			++m_nCurPos;

		if (m_nCurPos >= m_cs.GetLength())
			return FALSE;

		int nStartPos = m_nCurPos;
		while (m_nCurPos < m_cs.GetLength() && !m_sDelimiter[static_cast<BYTE>(m_cs[m_nCurPos])])
			++m_nCurPos;

		cs = m_cs.Mid(nStartPos, m_nCurPos - nStartPos);

		return TRUE;
	}

	CString Tail() const
	{
		int nCurPos = m_nCurPos;

		while (nCurPos < m_cs.GetLength() && m_sDelimiter[static_cast<BYTE>(m_cs[nCurPos])])
			++nCurPos;

		CString csResult;

		if (nCurPos < m_cs.GetLength())
			csResult = m_cs.Mid(nCurPos);

		return csResult;
	}

private:
	CString m_cs;
	std::bitset<256> m_sDelimiter;
	int m_nCurPos;
};
//
//--) // class CTokenizer
//
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
//
//(-- class CMatrix
//
class CMatrix
{
	//
	// 公有接口函数
	//
public:
	//
	// 构造与析构
	//

	CMatrix(); // 基础构造函数

	CMatrix(int nRows, int nCols); // 指定行列构造函数

	CMatrix(int nRows, int nCols, double value[]); // 指定数据构造函数

	CMatrix(int nSize); // 方阵构造函数

	CMatrix(int nSize, double value[]); // 指定数据方阵构造函数

	CMatrix(const CMatrix& other); // 拷贝构造函数

	BOOL Init(int nRows, int nCols); // 初始化矩阵

	BOOL MakeUnitMatrix(int nSize); // 将方阵初始化为单位矩阵

	virtual ~CMatrix(); // 析构函数

	//
	// 输入与显示
	//

	// 将字符串转换为矩阵数据
	BOOL FromString(CString s, const CString& sDelim = " ", BOOL bLineBreak = TRUE);
	// 将矩阵转换为字符串
	CString ToString(const CString& sDelim = " ", BOOL bLineBreak = TRUE) const;
	// 将矩阵的指定行转换为字符串
	CString RowToString(int nRow, const CString& sDelim = " ") const;
	// 将矩阵的指定列转换为字符串
	CString ColToString(int nCol, const CString& sDelim = " ") const;

	//
	// 元素与值操作
	//

	BOOL SetElement(int nRow, int nCol, double value); // 设置指定元素的值

	double GetElement(int nRow, int nCol) const; // 获取指定元素的值

	void SetData(double value[]); // 设置矩阵的值

	int GetNumColumns() const; // 获取矩阵的列数

	int GetNumRows() const; // 获取矩阵的行数

	int GetRowVector(int nRow, double* pVector) const; // 获取矩阵的指定行矩阵

	int GetColVector(int nCol, double* pVector) const; // 获取矩阵的指定列矩阵

	double* GetData() const; // 获取矩阵的值

	//
	// 数学操作
	//

	CMatrix& operator=(const CMatrix& other);
	BOOL operator==(const CMatrix& other) const;
	BOOL operator!=(const CMatrix& other) const;
	CMatrix operator+(const CMatrix& other) const;
	CMatrix operator-(const CMatrix& other) const;
	CMatrix operator*(double value) const;
	CMatrix operator*(const CMatrix& other) const;

	// 矩阵的转置
	CMatrix Transpose() const;

	//
	// 算法
	//

	// 实矩阵求逆的全选主元高斯－约当法
	BOOL InvertGaussJordan();

	// 对称正定矩阵的求逆
	BOOL InvertSsgj();

	//
	// 保护性数据成员
	//

protected:
	int m_nNumColumns; // 矩阵列数
	int m_nNumRows;	   // 矩阵行数
	double* m_pData;   // 矩阵数据缓冲区

	//
	// 内部函数
	//
private:
	void ppp(double a[], double e[], double s[], double v[], int m, int n);
	void sss(double fg[2], double cs[2]);
};
//
//--) // class CMatrix
//
//////////////////////////////////////////////////////////////////////////////////////
