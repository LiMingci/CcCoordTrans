//////////////////////////////////////////////////////////////////////
// Matrix.h
//
// ����������� CMatrix �������ӿ�
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
	// ���нӿں���
	//
public:
	//
	// ����������
	//

	CMatrix(); // �������캯��

	CMatrix(int nRows, int nCols); // ָ�����й��캯��

	CMatrix(int nRows, int nCols, double value[]); // ָ�����ݹ��캯��

	CMatrix(int nSize); // �����캯��

	CMatrix(int nSize, double value[]); // ָ�����ݷ����캯��

	CMatrix(const CMatrix& other); // �������캯��

	BOOL Init(int nRows, int nCols); // ��ʼ������

	BOOL MakeUnitMatrix(int nSize); // �������ʼ��Ϊ��λ����

	virtual ~CMatrix(); // ��������

	//
	// ��������ʾ
	//

	// ���ַ���ת��Ϊ��������
	BOOL FromString(CString s, const CString& sDelim = " ", BOOL bLineBreak = TRUE);
	// ������ת��Ϊ�ַ���
	CString ToString(const CString& sDelim = " ", BOOL bLineBreak = TRUE) const;
	// �������ָ����ת��Ϊ�ַ���
	CString RowToString(int nRow, const CString& sDelim = " ") const;
	// �������ָ����ת��Ϊ�ַ���
	CString ColToString(int nCol, const CString& sDelim = " ") const;

	//
	// Ԫ����ֵ����
	//

	BOOL SetElement(int nRow, int nCol, double value); // ����ָ��Ԫ�ص�ֵ

	double GetElement(int nRow, int nCol) const; // ��ȡָ��Ԫ�ص�ֵ

	void SetData(double value[]); // ���þ����ֵ

	int GetNumColumns() const; // ��ȡ���������

	int GetNumRows() const; // ��ȡ���������

	int GetRowVector(int nRow, double* pVector) const; // ��ȡ�����ָ���о���

	int GetColVector(int nCol, double* pVector) const; // ��ȡ�����ָ���о���

	double* GetData() const; // ��ȡ�����ֵ

	//
	// ��ѧ����
	//

	CMatrix& operator=(const CMatrix& other);
	BOOL operator==(const CMatrix& other) const;
	BOOL operator!=(const CMatrix& other) const;
	CMatrix operator+(const CMatrix& other) const;
	CMatrix operator-(const CMatrix& other) const;
	CMatrix operator*(double value) const;
	CMatrix operator*(const CMatrix& other) const;

	// �����ת��
	CMatrix Transpose() const;

	//
	// �㷨
	//

	// ʵ���������ȫѡ��Ԫ��˹��Լ����
	BOOL InvertGaussJordan();

	// �Գ��������������
	BOOL InvertSsgj();

	//
	// ���������ݳ�Ա
	//

protected:
	int m_nNumColumns; // ��������
	int m_nNumRows;	   // ��������
	double* m_pData;   // �������ݻ�����

	//
	// �ڲ�����
	//
private:
	void ppp(double a[], double e[], double s[], double v[], int m, int n);
	void sss(double fg[2], double cs[2]);
};
//
//--) // class CMatrix
//
//////////////////////////////////////////////////////////////////////////////////////
