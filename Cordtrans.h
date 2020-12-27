#pragma once
#include "Matrix.h"

struct PointNum //��������ݽṹ��
{
	char Pname[20]; //����
	double X;		//X����
	double Y;		//Y����
	double Z;		//Z����
};

//��������
//һ������㼯����Ե㼯�е����ز���
class CCrdNum
{
public:
	CCrdNum() //���캯����
	{
		memNum = 0;
	}
	CCrdNum(const CCrdNum& c) //���ƹ��캯��
	{
		memNum = c.memNum;				 //���Ƶ����
		Coord = new PointNum[memNum];	 //���Ƶ���
		for (int i = 0; i < memNum; i++) //���Ƶ�����
		{
			strcpy(Coord[i].Pname, c.Coord[i].Pname);
			Coord[i].X = c.Coord[i].X;
			Coord[i].Y = c.Coord[i].Y;
			Coord[i].Z = c.Coord[i].Z;
		}
	}

	CCrdNum& operator=(const CCrdNum& c)
	{
		if (this == &c)
		{
			return *this;
		}
		memNum = c.memNum;				 //���Ƶ����
		Coord = new PointNum[memNum];	 //���Ƶ���
		for (int i = 0; i < memNum; i++) //���Ƶ�����
		{
			strcpy(Coord[i].Pname, c.Coord[i].Pname);
			Coord[i].X = c.Coord[i].X;
			Coord[i].Y = c.Coord[i].Y;
			Coord[i].Z = c.Coord[i].Z;
		}
		return *this;
	}
	void Init(int num) //��ʼ���㼯
	{
		memNum = num;
		Coord = new PointNum[num];
	}

	~CCrdNum() //��������
	{
		if (Coord)
		{
			delete[] Coord;
			Coord = NULL;
		}
	}

	//��ȡ�㼯�е�ĸ���
	int GetmemNum()
	{
		return memNum;
	}

	//�����k����ĵ���������
	void Setdata(int k, char* name, char* X, char* Y, char* Z)
	{
		if (k > memNum)
		{
			// AfxMessageBox("�������ݴ���");
			return;
		}
		else
		{
			strcpy(Coord[k].Pname, name);
			Coord[k].X = atof(X);
			Coord[k].Y = atof(Y);
			Coord[k].Z = atof(Z);
		}
	}

	//�����k����ĵ���������
	void Setdata(int k, char* name, double X, double Y, double Z)
	{
		if (k > memNum)
		{
			// AfxMessageBox("�������ݴ���");
			return;
		}
		else
		{
			strcpy(Coord[k].Pname, name);
			Coord[k].X = X;
			Coord[k].Y = Y;
			Coord[k].Z = Z;
		}
	}

	//��ȡ��k����ĵ���
	char* Getdata(int k)
	{
		return Coord[k].Pname;
	}

	//��ȡ��k�����ĳ��������ֵ
	//index=1Ϊ��ȡX����
	//index=2Ϊ��ȡY����
	//index=3Ϊ��ȡZ����
	double Getdata(int k, int index)
	{
		double data;
		switch (index)
		{
		case 1:
			data = Coord[k].X;
			break;
		case 2:
			data = Coord[k].Y;
			break;
		case 3:
			data = Coord[k].Z;
			break;
		}
		return data;
	}

	//��������֮��ľ���
	//index1��index2�ֱ�Ϊ������
	double Distance(int index1, int index2)
	{
		double distance = 0.0;
		if ((index1 < memNum) && (index2 < memNum))
		{

			distance = sqrt(pow((Coord[index1].X - Coord[index2].X), 2) + pow((Coord[index1].Y - Coord[index2].Y), 2) + pow((Coord[index1].Z - Coord[index2].Z), 2));
		}
		else
		{
			// AfxMessageBox("�±�Խ��",MB_OK|MB_ICONINFORMATION);
			distance = -1;
		}

		return distance;
	}

	//���ַ�������ʽ����㼯�����е����Ϣ
	CString Printpoint()
	{
		CString txt = "\n";
		CString crdX, crdY, crdZ;
		txt = txt + "����" + "          " + "X" + "          " + "Y" + "           " + "Z" + "\n";
		for (int i = 0; i < memNum; i++)
		{
			crdX.Format("%0.4lf", Getdata(i, 1));
			crdY.Format("%0.4lf", Getdata(i, 2));
			crdZ.Format("%0.4lf", Getdata(i, 3));
			txt = txt + Getdata(i) + "     " + crdX + "     " + crdY + "     " + crdZ + "\n";
		}
		return txt;
	}

private:
	int memNum;		 //�����
	PointNum* Coord; //������
};

//����ת����
//ʵ�ּ�����ת����ƽ�Ʋ���
//ʵ��˳��ƥ�䣬����ƥ�䣬����ƥ�䣬����ƥ����㷽ʽ
//ʵ�������㼯֮��Ѱ��ͬ����
//ʵ�������㼯�������ȡ���
class CCordtrans
{
public:
	//���캯��
	CCordtrans(CCrdNum num1, CCrdNum num2);
	CCordtrans();
	//��������
	virtual ~CCordtrans();

	CString FreeMatch(CMatrix& R, CMatrix& T, CCrdNum& Coord1to2);
	CString InitialMatch(CMatrix& R, CMatrix& T, CCrdNum& Coord1to2);
	CString PointNameMatch(CMatrix& R, CMatrix& T, CCrdNum& Coord1to2);
	CString DistanceMatch(CMatrix& R, CMatrix& T, CCrdNum& Coord1to2);
	//������ת��
	CString AngleofSwing(CMatrix R0);
	//���Ϊ�ַ���
	CString Print(CMatrix R0, CMatrix T0, CCrdNum crd3, CCrdNum crd4, CCrdNum crd1to2);
	//��ȷ�������
	void PrecisionCalPar(CMatrix& RX, CMatrix& TX);
	//��ȷ����ת������
	void PrecisionCalPar(CCrdNum crd1, CCrdNum crd2, CMatrix& RX, CMatrix& TX);
	//���㷽���
	CMatrix DirectionAngle(CMatrix R0);
	//�����ļ�
	void Save();
	//�����λƫ��
	CString Differencebetween(CCrdNum crd1, CCrdNum crd2);
	//����ƥ��
	CString FreeMatch(CCrdNum& Coord1to2);
	//��ȡ�ĸ������
	void DistillPoint(CCrdNum& crd3, CCrdNum& crd4);
	//Ѱ�Ҿ������һ����ֵ�ĵ��
	void FindDistancePoint(int FIGURE, CCrdNum& crd3, CCrdNum& crd4);
	//������ֵ
	double CalThreshold(CCrdNum crd1, CCrdNum crd2);
	//����ƥ��
	CString DistanceMatch(CCrdNum& Coord1to2);
	//˳��ƥ��
	CString InitialMatch(CCrdNum& Coord1to2);
	//����ƥ��
	CString PointNameMatch(CCrdNum& Coord1to2);
	//Ѱ��ͬ�����
	void FindSameName(CCrdNum crd1, CCrdNum crd2, CCrdNum& crd3, CCrdNum& crd4);
	//Ѱ��ͬ�����
	void FindSameName(CCrdNum& crd1, CCrdNum& crd2);
	//����ת������
	void CalPar(CCrdNum crd1, CCrdNum crd2, CMatrix& RX, CMatrix& TX);
	void CalPar(CMatrix& RX, CMatrix& TX);
	//����ת������
	void Translate(CCrdNum& Coord1to2);
	void Translate(CMatrix R0, CMatrix T0, CCrdNum crd1, CCrdNum& Coord1to2);

public:
	//������ת��
	static void AngleofSwing(CMatrix R0, CString& Rx, CString& Ry, CString& Rz);
	static CString AngleofSwing(CMatrix R0, double Rx, double Ry, double Rz);
	//����ת��Ϊ��
	static double R2D(double radian);

private:
	int Ncount1;	//���������
	int Ncount2;	//�ο������
	CCrdNum Coord1; //������㼯
	CCrdNum Coord2; //�ο���㼯
	CMatrix R;		//��ת����
	CMatrix T;		//ƽ�Ʋ���
	CString txt;	//����ַ�������
};
