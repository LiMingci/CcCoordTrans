#include "����ת��.h"
#include "Openfile.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//���캯��
//pathΪ�ļ�·����captionnameΪ��¼�ļ��Ƿ��б���
//pointnameΪ��¼�ļ����Ƿ��е�����divnameΪ�ָ������
COpenfile::COpenfile(CString path,int captionname,int pointname,int divname)
{
	cname=captionname;
	poname=pointname;
	m_pathname=path;
	divmark=divname;
	filedata="";
	Ncount=0;
	File.Open(m_pathname,CFile::modeRead);//��ֻ����ʽ���ļ�

}

//��������
COpenfile::~COpenfile()
{
	File.Close();//�ر��ļ�
	Ncount=0;//������0

}


//��ȡ�ļ����ж���������
int COpenfile::Findncount()
{
	int N=0;
	DWORD dwActual=File.SeekToEnd();//����ļ�ĩβ��ָ��
	File.SeekToBegin();//���ļ�ָ������ָ���ļ�ͷ
	if (cname)//�ж��Ƿ��б���ͷ
	{
		File.ReadString(filedata);
	}
	
	do //�����Ч�ļ�����
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
	File.SeekToBegin();//��ָ������ָ���ļ�ͷ
	return N;
}

//��ȡ�ļ���k����Ϣ
//nameΪ������tempXΪX���꣬tempYΪY���꣬tempZΪZ����
void COpenfile::Readdata(int k,char *name, char *tempX, char *tempY, char *tempZ)
{
	
	CString number;//�޵���ʱ�����������
  	if (cname)//�ж����ޱ���ͷ
  	{
		File.ReadString(filedata);
		cname=!cname;
  	}
	
	do //�ж��Ƿ��ǿ��в�ɾ������
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
	
	switch (divmark)//�жϷָ�����Ͳ�����Ӧ��ȡ
	{
	case 1:
		if (poname==0)//���û�е���
		{
			number.Format("%d",k);
			number="Point"+number;
			strcpy(name,number);
			sscanf(filedata,"%s %s %s",tempX,tempY,tempZ);
			
		}
		else //����е���
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
