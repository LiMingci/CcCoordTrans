#pragma once
#include "Matrix.h"

struct PointNum //定义点数据结构体
{
	char Pname[20]; //点名
	double X;		//X坐标
	double Y;		//Y坐标
	double Z;		//Z坐标
};

//坐标组类
//一个坐标点集及其对点集中点的相关操作
class CCrdNum
{
public:
	CCrdNum() //构造函数对
	{
		memNum = 0;
	}
	CCrdNum(const CCrdNum& c) //复制构造函数
	{
		memNum = c.memNum;				 //复制点个数
		Coord = new PointNum[memNum];	 //复制点名
		for (int i = 0; i < memNum; i++) //复制点坐标
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
		memNum = c.memNum;				 //复制点个数
		Coord = new PointNum[memNum];	 //复制点名
		for (int i = 0; i < memNum; i++) //复制点坐标
		{
			strcpy(Coord[i].Pname, c.Coord[i].Pname);
			Coord[i].X = c.Coord[i].X;
			Coord[i].Y = c.Coord[i].Y;
			Coord[i].Z = c.Coord[i].Z;
		}
		return *this;
	}
	void Init(int num) //初始化点集
	{
		memNum = num;
		Coord = new PointNum[num];
	}

	~CCrdNum() //析构函数
	{
		if (Coord)
		{
			delete[] Coord;
			Coord = NULL;
		}
	}

	//获取点集中点的个数
	int GetmemNum()
	{
		return memNum;
	}

	//输入第k个点的点名和坐标
	void Setdata(int k, char* name, char* X, char* Y, char* Z)
	{
		if (k > memNum)
		{
			// AfxMessageBox("输入数据错误");
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

	//输入第k个点的点名和坐标
	void Setdata(int k, char* name, double X, double Y, double Z)
	{
		if (k > memNum)
		{
			// AfxMessageBox("输入数据错误");
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

	//获取第k个点的点名
	char* Getdata(int k)
	{
		return Coord[k].Pname;
	}

	//获取第k个点的某个轴坐标值
	//index=1为获取X坐标
	//index=2为获取Y坐标
	//index=3为获取Z坐标
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

	//计算两点之间的距离
	//index1，index2分别为两点点号
	double Distance(int index1, int index2)
	{
		double distance = 0.0;
		if ((index1 < memNum) && (index2 < memNum))
		{

			distance = sqrt(pow((Coord[index1].X - Coord[index2].X), 2) + pow((Coord[index1].Y - Coord[index2].Y), 2) + pow((Coord[index1].Z - Coord[index2].Z), 2));
		}
		else
		{
			// AfxMessageBox("下标越界",MB_OK|MB_ICONINFORMATION);
			distance = -1;
		}

		return distance;
	}

	//以字符串的形式输出点集中所有点的信息
	CString Printpoint()
	{
		CString txt = "\n";
		CString crdX, crdY, crdZ;
		txt = txt + "点名" + "          " + "X" + "          " + "Y" + "           " + "Z" + "\n";
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
	int memNum;		 //点个数
	PointNum* Coord; //点数组
};

//坐标转化类
//实现计算旋转参数平移参数
//实现顺序匹配，点名匹配，自由匹配，距离匹配计算方式
//实现在两点集之间寻找同名点
//实现在两点集中随机抽取点对
class CCordtrans
{
public:
	//构造函数
	CCordtrans(CCrdNum num1, CCrdNum num2);
	CCordtrans();
	//析构函数
	virtual ~CCordtrans();

	CString FreeMatch(CMatrix& R, CMatrix& T, CCrdNum& Coord1to2);
	CString InitialMatch(CMatrix& R, CMatrix& T, CCrdNum& Coord1to2);
	CString PointNameMatch(CMatrix& R, CMatrix& T, CCrdNum& Coord1to2);
	CString DistanceMatch(CMatrix& R, CMatrix& T, CCrdNum& Coord1to2);
	//计算旋转角
	CString AngleofSwing(CMatrix R0);
	//输出为字符串
	CString Print(CMatrix R0, CMatrix T0, CCrdNum crd3, CCrdNum crd4, CCrdNum crd1to2);
	//精确计算参数
	void PrecisionCalPar(CMatrix& RX, CMatrix& TX);
	//精确计算转换参数
	void PrecisionCalPar(CCrdNum crd1, CCrdNum crd2, CMatrix& RX, CMatrix& TX);
	//计算方向角
	CMatrix DirectionAngle(CMatrix R0);
	//保存文件
	void Save();
	//计算点位偏差
	CString Differencebetween(CCrdNum crd1, CCrdNum crd2);
	//自由匹配
	CString FreeMatch(CCrdNum& Coord1to2);
	//提取四个随机点
	void DistillPoint(CCrdNum& crd3, CCrdNum& crd4);
	//寻找距离符合一定阈值的点对
	void FindDistancePoint(int FIGURE, CCrdNum& crd3, CCrdNum& crd4);
	//计算阈值
	double CalThreshold(CCrdNum crd1, CCrdNum crd2);
	//距离匹配
	CString DistanceMatch(CCrdNum& Coord1to2);
	//顺序匹配
	CString InitialMatch(CCrdNum& Coord1to2);
	//点名匹配
	CString PointNameMatch(CCrdNum& Coord1to2);
	//寻找同名点对
	void FindSameName(CCrdNum crd1, CCrdNum crd2, CCrdNum& crd3, CCrdNum& crd4);
	//寻找同名点对
	void FindSameName(CCrdNum& crd1, CCrdNum& crd2);
	//计算转换参数
	void CalPar(CCrdNum crd1, CCrdNum crd2, CMatrix& RX, CMatrix& TX);
	void CalPar(CMatrix& RX, CMatrix& TX);
	//计算转换坐标
	void Translate(CCrdNum& Coord1to2);
	void Translate(CMatrix R0, CMatrix T0, CCrdNum crd1, CCrdNum& Coord1to2);

public:
	//计算旋转角
	static void AngleofSwing(CMatrix R0, CString& Rx, CString& Ry, CString& Rz);
	static CString AngleofSwing(CMatrix R0, double Rx, double Ry, double Rz);
	//弧度转化为度
	static double R2D(double radian);

private:
	int Ncount1;	//测量点个数
	int Ncount2;	//参考点个数
	CCrdNum Coord1; //测量点点集
	CCrdNum Coord2; //参考点点集
	CMatrix R;		//旋转矩阵
	CMatrix T;		//平移参数
	CString txt;	//输出字符串变量
};
