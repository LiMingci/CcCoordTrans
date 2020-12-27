#pragma once

//打开文件并获取文件信息
class COpenfile
{
public:
	void Readdata(int k, char* name, char* tempX, char* tempY, char* tempZ); //读取文件第k行数据
	int Findncount();														 //计算文件行数
	COpenfile(CString path, int captionname, int pointname, int divname);	 //构造函数
	virtual ~COpenfile();													 //析构函数
private:
	CString m_pathname; //文件路径
	int divmark;		//数据之间的分隔符
	int cname;			//判断有无标题
	int poname;			//判断有无点名
	int Ncount;			//文件行数
	CString filedata;	//临时存取变量
	CStdioFile File;
};
