#pragma once

//���ļ�����ȡ�ļ���Ϣ
class COpenfile
{
public:
	void Readdata(int k, char* name, char* tempX, char* tempY, char* tempZ); //��ȡ�ļ���k������
	int Findncount();														 //�����ļ�����
	COpenfile(CString path, int captionname, int pointname, int divname);	 //���캯��
	virtual ~COpenfile();													 //��������
private:
	CString m_pathname; //�ļ�·��
	int divmark;		//����֮��ķָ���
	int cname;			//�ж����ޱ���
	int poname;			//�ж����޵���
	int Ncount;			//�ļ�����
	CString filedata;	//��ʱ��ȡ����
	CStdioFile File;
};
