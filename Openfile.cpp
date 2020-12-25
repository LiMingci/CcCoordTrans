#include "坐标转换.h"
#include "Openfile.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//构造函数
//path为文件路径，captionname为记录文件是否有标题
//pointname为记录文件中是否有点名，divname为分割符类型
COpenfile::COpenfile(CString path,int captionname,int pointname,int divname)
{
	cname=captionname;
	poname=pointname;
	m_pathname=path;
	divmark=divname;
	filedata="";
	Ncount=0;
	File.Open(m_pathname,CFile::modeRead);//以只读方式打开文件

}

//析构函数
COpenfile::~COpenfile()
{
	File.Close();//关闭文件
	Ncount=0;//行数置0

}


//获取文件中有多少行数据
int COpenfile::Findncount()
{
	int N=0;
	DWORD dwActual=File.SeekToEnd();//获得文件末尾的指针
	File.SeekToBegin();//将文件指针重新指向文件头
	if (cname)//判断是否有标题头
	{
		File.ReadString(filedata);
	}
	
	do //获得有效文件行数
	{
		
		if(File.ReadString(filedata))
		{
			filedata.TrimLeft();
			filedata.TrimRight();
		}
		if (filedata=="")
		{
			Ncount=Ncount-1;
			
		}
		Ncount=Ncount+1;
	} while (File.GetPosition()!=dwActual);
	N=Ncount;
	File.SeekToBegin();//将指针重新指向文件头
	return N;
}

//读取文件第k行信息
//name为点名，tempX为X坐标，tempY为Y坐标，tempZ为Z坐标
void COpenfile::Readdata(int k,char *name, char *tempX, char *tempY, char *tempZ)
{
	
	CString number;//无点名时替代点名变量
  	if (cname)//判断有无标题头
  	{
		File.ReadString(filedata);
		cname=!cname;
  	}
	
	do //判断是否是空行并删除空行
	{
		if(File.ReadString(filedata))
		{
			filedata.TrimLeft();
			filedata.TrimRight();
			
		}
		if (filedata=="")
		{
			//File.Seek(0,CFile::current);
			File.ReadString(filedata);
			
		}
	} while (filedata.IsEmpty());
	
	switch (divmark)//判断分割符类型并作相应读取
	{
	case 1:
		if (poname==0)//如果没有点名
		{
			number.Format("%d",k);
			number="Point"+number;
			strcpy(name,number);
			sscanf(filedata,"%s %s %s",tempX,tempY,tempZ);
			
		}
		else //如果有点名
			sscanf(filedata,"%s %s %s %s",name,tempX,tempY,tempZ);
		break;
	case 2:
		if (poname==0)
		{
			number.Format("%d",k);
			number="Point"+number;
			strcpy(name,number);
			sscanf(filedata,"%[^','],%[^','],%s",tempX,tempY,tempZ);
			
		}
		else
			sscanf(filedata,"%[^','],%[^','],%[^','],%s",name,tempX,tempY,tempZ);
		break;
	case 3:
		if (poname==0)
		{
			number.Format("%d",k);
			number="Point"+number;
			strcpy(name,number);
			sscanf(filedata,"%s %s %s",tempX,tempY,tempZ);
			
		}
		else
			sscanf(filedata,"%s %s %s %s",name,tempX,tempY,tempZ);
		break;
	case 4:
		if (poname==0)
		{
			number.Format("%d",k);
			number="Point"+number;
			strcpy(name,number);
			sscanf(filedata,"%[^';'];%[^';'];%s",tempX,tempY,tempZ);
			
		}
		else
			sscanf(filedata,"%[^';'];%[^';'];%[^';'];%s",name,tempX,tempY,tempZ);
		break;
	case 5:
		if (poname==0)
		{
			number.Format("%d",k);
			number="Point"+number;
			strcpy(name,number);
			sscanf(filedata,"%[^':']:%[^':']:%s",tempX,tempY,tempZ);
			
		}
		else
			sscanf(filedata,"%[^':']:%[^':']:%[^':']:%s",name,tempX,tempY,tempZ);
		break;
	}
	

}
