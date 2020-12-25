
#include "����ת��.h"
#include "Cordtrans.h"
#include "Openfile.h"
#include <string.h>
#define MEASURE 1
#define REFERENCE 2
#define PI 3.14159265359

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//���ƹ��캯��
//��ʼ��R��ת�����Tƽ�Ʋ���
//��ʼ��������㼯Coord1
//��ʼ���ο���㼯Coord2
 CCordtrans::CCordtrans(CCrdNum num1,CCrdNum num2):R(3,3),T(3,1),Coord1(num1),Coord2(num2)//���캯����ʼ��
 {
 	
 	Ncount1=num1.GetmemNum();
 	Ncount2=num2.GetmemNum();
 }

 //���캯��
 CCordtrans::CCordtrans():R(3,3),T(3,1)
 {
	txt="\n";
	Ncount1=Ncount2=0;

 }
//��������
CCordtrans::~CCordtrans()
{
	Ncount1=0;
	Ncount1=0;
	txt="\n";


}

////����ת������
void CCordtrans::CalPar(CMatrix &RX, CMatrix &TX)
{
	
	CMatrix M(3*Ncount2,6);//ϵ������
	CMatrix D(6,1);//��ز���
	CMatrix L(3*Ncount2,1);//������
	CMatrix InvM(6,6);
	double a,b,c;
	//���ϵ������ͳ�������
	for(int i=0;i<Ncount2;i++)
	{
		
		M.SetElement(3*i,0,0.00);
		M.SetElement(3*i,1,(-Coord2.Getdata(i,3)-Coord1.Getdata(i,3)));
		M.SetElement(3*i,2,(Coord2.Getdata(i,2)+Coord1.Getdata(i,2)));
		M.SetElement(3*i,3,1.00);
		M.SetElement(3*i,4,0.00);
		M.SetElement(3*i,5,0.00);
		M.SetElement(3*i+1,0,(Coord2.Getdata(i,3)+Coord1.Getdata(i,3)));
		M.SetElement(3*i+1,1,0.00);
		M.SetElement(3*i+1,2,(-Coord2.Getdata(i,1)-Coord1.Getdata(i,1)));
		M.SetElement(3*i+1,3,0.00);
		M.SetElement(3*i+1,4,1.00);
		M.SetElement(3*i+1,5,0.00);
		M.SetElement(3*i+2,0,(-Coord2.Getdata(i,2)-Coord1.Getdata(i,2)));
		M.SetElement(3*i+2,1,(Coord2.Getdata(i,1)+Coord1.Getdata(i,1)));
		M.SetElement(3*i+2,2,0.00);
		M.SetElement(3*i+2,3,0.00);
		M.SetElement(3*i+2,4,0.00);
		M.SetElement(3*i+2,5,1.00);
		L.SetElement(3*i,0,(Coord1.Getdata(i,1)-Coord2.Getdata(i,1)));
		L.SetElement(3*i+1,0,(Coord1.Getdata(i,2)-Coord2.Getdata(i,2)));
		L.SetElement(3*i+2,0,(Coord1.Getdata(i,3)-Coord2.Getdata(i,3)));
	}
	InvM=M.Transpose()*M;
	if(InvM.InvertGaussJordan())
	{
		D=InvM*(M.Transpose()*L);
	}

	a=D.GetElement(0,0);
	b=D.GetElement(1,0);
	c=D.GetElement(2,0);
	CMatrix I_S(3,3);//���������м����
	CMatrix U(3,1);//���������м����
	I_S.SetElement(0,0,1);
	I_S.SetElement(0,1,c);
	I_S.SetElement(0,2,-b);
	I_S.SetElement(1,0,-c);
	I_S.SetElement(1,1,1);
	I_S.SetElement(1,2,a);
	I_S.SetElement(2,0,b);
	I_S.SetElement(2,1,-a);
	I_S.SetElement(2,2,1);
	U.SetElement(0,0,D.GetElement(3,0));
	U.SetElement(1,0,D.GetElement(4,0));
	U.SetElement(2,0,D.GetElement(5,0));
	if(I_S.InvertGaussJordan())
	{
		T=(I_S*U)*(-1);//���ƽ�Ʋ���
	}
	//������ת����
	R.SetElement(0,0,(1+a*a-b*b-c*c));
	R.SetElement(0,1,(2*(a*b-c)));
	R.SetElement(0,2,(2*(a*c+b)));
	R.SetElement(1,0,(2*(a*b+c)));
	R.SetElement(1,1,(1-a*a+b*b-c*c));
	R.SetElement(1,2,(2*(b*c-a)));
	R.SetElement(2,0,(2*(a*c-b)));
	R.SetElement(2,1,(2*(b*c+a)));
	R.SetElement(2,2,(1-a*a-b*b+c*c));
 	R=R*(1/(1+a*a+b*b+c*c));

	RX=R;
	TX=T;
	
}

//����ת������

void CCordtrans::CalPar(CCrdNum crd1, CCrdNum crd2, CMatrix &RX, CMatrix &TX)
{
	CMatrix M(3*crd2.GetmemNum(),6);//ϵ������
	CMatrix D(6,1);//��ز���
	CMatrix L(3*crd2.GetmemNum(),1);//������
	CMatrix InvM(6,6);
	double a,b,c;
	for(int i=0;i<crd2.GetmemNum();i++)//���ϵ������ͳ�������
	{
		
		M.SetElement(3*i,0,0.00);
		M.SetElement(3*i,1,(-crd2.Getdata(i,3)-crd1.Getdata(i,3)));
		M.SetElement(3*i,2,(crd2.Getdata(i,2)+crd1.Getdata(i,2)));
		M.SetElement(3*i,3,1.00);
		M.SetElement(3*i,4,0.00);
		M.SetElement(3*i,5,0.00);
		M.SetElement(3*i+1,0,(crd2.Getdata(i,3)+crd1.Getdata(i,3)));
		M.SetElement(3*i+1,1,0.00);
		M.SetElement(3*i+1,2,(-crd2.Getdata(i,1)-crd1.Getdata(i,1)));
		M.SetElement(3*i+1,3,0.00);
		M.SetElement(3*i+1,4,1.00);
		M.SetElement(3*i+1,5,0.00);
		M.SetElement(3*i+2,0,(-crd2.Getdata(i,2)-crd1.Getdata(i,2)));
		M.SetElement(3*i+2,1,(crd2.Getdata(i,1)+crd1.Getdata(i,1)));
		M.SetElement(3*i+2,2,0.00);
		M.SetElement(3*i+2,3,0.00);
		M.SetElement(3*i+2,4,0.00);
		M.SetElement(3*i+2,5,1.00);
		L.SetElement(3*i,0,(crd1.Getdata(i,1)-crd2.Getdata(i,1)));
		L.SetElement(3*i+1,0,(crd1.Getdata(i,2)-crd2.Getdata(i,2)));
		L.SetElement(3*i+2,0,(crd1.Getdata(i,3)-crd2.Getdata(i,3)));
	}
	InvM=M.Transpose()*M;
	if(InvM.InvertGaussJordan())
	{
		D=InvM*(M.Transpose()*L);
	}
	
	a=D.GetElement(0,0);
	b=D.GetElement(1,0);
	c=D.GetElement(2,0);
	CMatrix I_S(3,3);
	CMatrix U(3,1);//���������м����
	I_S.SetElement(0,0,1);
	I_S.SetElement(0,1,c);
	I_S.SetElement(0,2,-b);
	I_S.SetElement(1,0,-c);
	I_S.SetElement(1,1,1);
	I_S.SetElement(1,2,a);
	I_S.SetElement(2,0,b);
	I_S.SetElement(2,1,-a);
	I_S.SetElement(2,2,1);
	U.SetElement(0,0,D.GetElement(3,0));
	U.SetElement(1,0,D.GetElement(4,0));
	U.SetElement(2,0,D.GetElement(5,0));
	if(I_S.InvertGaussJordan())
	{
		T=(I_S*U)*(-1);//���ƽ�Ʋ���
	}
	R.SetElement(0,0,(1+a*a-b*b-c*c));//������ת����
	R.SetElement(0,1,(2*(a*b-c)));
	R.SetElement(0,2,(2*(a*c+b)));
	R.SetElement(1,0,(2*(a*b+c)));
	R.SetElement(1,1,(1-a*a+b*b-c*c));
	R.SetElement(1,2,(2*(b*c-a)));
	R.SetElement(2,0,(2*(a*c-b)));
	R.SetElement(2,1,(2*(b*c+a)));
	R.SetElement(2,2,(1-a*a-b*b+c*c));
	R=R*(1/(1+a*a+b*b+c*c));
	
	RX=R;
	TX=T;

	
}


void CCordtrans::Translate(CCrdNum &Coord1to2)//����ת������
{
	Coord1to2.Init(Ncount1);
	CMatrix CoordM1(3,Ncount1);//�������������
	CMatrix CoordM2(3,Ncount1);//������ת�����������
	CMatrix TM(3,Ncount1);//ת���м����
	double T1,T2,T3;//����ƽ�Ʋ���
	T1=T.GetElement(0,0);//��ƽ�Ʋ���������ȡֵ
	T2=T.GetElement(1,0);
	T3=T.GetElement(2,0);
	for(int i=0;i<Ncount1;i++)//���������������ֵ
	{
		CoordM1.SetElement(0,i,Coord1.Getdata(i,1));
		CoordM1.SetElement(1,i,Coord1.Getdata(i,2));
		CoordM1.SetElement(2,i,Coord1.Getdata(i,3));
		TM.SetElement(0,i,T1);
		TM.SetElement(1,i,T2);
		TM.SetElement(2,i,T3);
	}
	CoordM2=R*CoordM1+TM;//����ת��
	for (int i=0;i<Ncount1;i++)
	{
		Coord1to2.Setdata(i,Coord1.Getdata(i),CoordM2.GetElement(0,i),CoordM2.GetElement(1,i),CoordM2.GetElement(2,i));

	}
	
}

void CCordtrans::Translate(CMatrix R0,CMatrix T0,CCrdNum crd1,CCrdNum &Coord1to2)//����ת������
{
	int count1=crd1.GetmemNum();
	Coord1to2.Init(count1);
	CMatrix CoordM1(3,count1);//�������������
	CMatrix CoordM2(3,count1);//������ת�����������
	CMatrix TM(3,count1);//ת���м����
	double T1,T2,T3;//����ƽ�Ʋ���
	T1=T0.GetElement(0,0);//��ƽ�Ʋ���������ȡֵ
	T2=T0.GetElement(1,0);
	T3=T0.GetElement(2,0);
	for(int i=0;i<count1;i++)//���������������ֵ
	{
		CoordM1.SetElement(0,i,crd1.Getdata(i,1));
		CoordM1.SetElement(1,i,crd1.Getdata(i,2));
		CoordM1.SetElement(2,i,crd1.Getdata(i,3));
		TM.SetElement(0,i,T1);
		TM.SetElement(1,i,T2);
		TM.SetElement(2,i,T3);
	}
	CoordM2=R0*CoordM1+TM;//����ת��
	for (int i=0;i<count1;i++)
	{
		Coord1to2.Setdata(i,crd1.Getdata(i),CoordM2.GetElement(0,i),CoordM2.GetElement(1,i),CoordM2.GetElement(2,i));
		
	}
}

//��ȷ����ת������

void CCordtrans::PrecisionCalPar(CCrdNum crd1, CCrdNum crd2, CMatrix &RX, CMatrix &TX)
{
	CMatrix R1(3,3),T1(3,1);
	CMatrix R2(3,3),T2(3,1);
	CCrdNum crd1to2;
	
	CalPar(crd1,crd2,R1,T1);
	Translate(R1,T1,crd1,crd1to2);
	CalPar(crd1to2,crd2,R2,T2);
	
	RX=R2*R1;
	TX=R2*T1+T2;
	
}


//��ȷ����ת������

void CCordtrans::PrecisionCalPar(CMatrix &RX, CMatrix &TX)
{
	CMatrix R1(3,3),T1(3,1);
	CMatrix R2(3,3),T2(3,1);
	CCrdNum Coord1to2;
	
	CalPar(Coord1,Coord2,R1,T1);
	Translate(R1,T1,Coord1,Coord1to2);
	CalPar(Coord1to2,Coord2,R2,T2);
	
	RX=R2*R1;
	TX=R2*T1+T2;
	
}

//Ѱ��ͬ����
void CCordtrans::Findsamename(CCrdNum &crd3, CCrdNum &crd4)
{
	int Pcount=0;
    CString tes;
	for(int i=0;i<Ncount2;i++)
	{
		for(int j=0;j<Ncount1;j++)
		{
			if (strcmp(Coord2.Getdata(i),Coord1.Getdata(j))==0)
			{
				Pcount=Pcount+1;
			}
		}
		
	}
	tes.Format("%d",Pcount);
	tes=tes+"����ƥ���ϵĵ���";
	AfxMessageBox(tes,MB_OK|MB_ICONINFORMATION);
	crd3.Init(Pcount);
	crd4.Init(Pcount);
	int m=0;
	for(int i=0;i<Ncount2;i++)
	{
		for(int j=0;j<Ncount1;j++)
		{
			if (strcmp(Coord2.Getdata(i),Coord1.Getdata(j))==0)
			{
				crd3.Setdata(m,Coord1.Getdata(j),Coord1.Getdata(j,1),Coord1.Getdata(j,2),Coord1.Getdata(j,3));
				crd4.Setdata(m,Coord2.Getdata(i),Coord2.Getdata(i,1),Coord2.Getdata(i,2),Coord2.Getdata(i,3));
				m=m+1;
			}
		}	
	}

}

void CCordtrans::Findsamename(CCrdNum crd1, CCrdNum crd2, CCrdNum &crd3, CCrdNum &crd4)//Ѱ��ͬ����
{
	int Pcount=0;
    CString tes;
	int count1=crd1.GetmemNum();
	int count2=crd2.GetmemNum();
	for(int i=0;i<count2;i++)
	{
		for(int j=0;j<count1;j++)
		{
			if (strcmp(crd2.Getdata(i),crd1.Getdata(j))==0)
			{
				Pcount=Pcount+1;
			}
		}
		
	}
	tes.Format("%d",Pcount);
	tes=tes+"����ƥ���ϵĵ���";
	AfxMessageBox(tes,MB_OK|MB_ICONINFORMATION);
	crd3.Init(Pcount);
	crd4.Init(Pcount);
	int m=0;
	for(int i=0;i<count2;i++)
	{
		for(int j=0;j<count1;j++)
		{
			if (strcmp(crd2.Getdata(i),crd1.Getdata(j))==0)
			{
				crd3.Setdata(m,crd1.Getdata(j),crd1.Getdata(j,1),crd1.Getdata(j,2),crd1.Getdata(j,3));
				crd4.Setdata(m,crd2.Getdata(i),crd2.Getdata(i,1),crd2.Getdata(i,2),crd2.Getdata(i,3));
				m=m+1;
			}
		}	
	}

}

void CCordtrans::Distillpoint(CCrdNum &crd3, CCrdNum &crd4)//��ȡ�ĸ������
{
	double edge1=0.0,edge2=0.0,edge3=0.0,edge4=0.0,edge5=0.0,edge6=0.0;
	double distance1=0.0,distance2=0.0,distance3=0.0,distance4=0.0,distance5=0.0,distance6=0.0;
	double dem=0.33;//???????????????????????????????
	CString tes;
	CCrdNum crd40;
	int M=4,n=-1;
	crd40.Init(Ncount2);
	crd4.Init(M);
	crd3.Init(M);
	for (int i=0;i<Ncount2;i++)
	{
		crd40.Setdata(i,"\0",0.0,0.0,0.0);

	}
	for (int i=0;i<M;i++)
	{
		crd3.Setdata(i,"\0",0.0,0.0,0.0);
		crd4.Setdata(i,"\0",0.0,0.0,0.0);

	}

	do
	{
		srand((unsigned)time(NULL));//�������������
		for (int i=0;i<M;i++)//�Ӳο����л���ĸ�����Ĳ�ͬԪ��
		{
			n=rand()%Ncount2;//����һ�����
			if(strcmp(crd40.Getdata(n),"\0")==0)//ȡ��
			{
				crd40.Setdata(n,"NOName",0.0,0.0,0.0);
				crd4.Setdata(i,Coord2.Getdata(n),Coord2.Getdata(n,1),Coord2.Getdata(n,2),Coord2.Getdata(n,3));

				//AfxMessageBox(crd4.Getdata(i) ,MB_OK|MB_ICONINFORMATION);
				
			}
			else
				i--;
			
		}
		edge1=crd4.Distance(0,1);
		edge2=crd4.Distance(0,2);
		edge3=crd4.Distance(1,2);
	}while(fabs(edge1-edge2)<1 && fabs(edge1-edge3)<1 && fabs(edge2-edge3)<1);
	edge4=crd4.Distance(3,0);
	edge5=crd4.Distance(3,1);
	edge6=crd4.Distance(3,2);
	for (int i=0;i<Ncount1;i++)
	{
		crd3.Setdata(0,Coord1.Getdata(i),Coord1.Getdata(i,1),Coord1.Getdata(i,2),Coord1.Getdata(i,3));
		for (int j=0;j<Ncount1;j++)
		{
			if (i==j)
			{
				continue;
			}
			else
			{
				distance1=Coord1.Distance(i,j);
				if ((fabs(distance1-edge1)<dem))
				{
					crd3.Setdata(1,Coord1.Getdata(j),Coord1.Getdata(j,1),Coord1.Getdata(j,2),Coord1.Getdata(j,3));
				}
				else
				{
					distance2=Coord1.Distance(i,j);
					if ((fabs(distance2-edge2)<dem))
					{
						crd3.Setdata(2,Coord1.Getdata(j),Coord1.Getdata(j,1),Coord1.Getdata(j,2),Coord1.Getdata(j,3));

					}

				}

			}

		}
		if(fabs(crd3.Getdata(1,1))>0 && fabs(crd3.Getdata(2,1))>0)
		{
			distance3=crd3.Distance(1,2);
			if (fabs(distance3-edge3)<dem)
			{
				break;
			}
			
		}

	}
	for (int i=0;i<Ncount1;i++)
	{
		distance4=sqrt(pow((Coord1.Getdata(i,1)-crd3.Getdata(0,1)),2)+pow((Coord1.Getdata(i,2)-crd3.Getdata(0,2)),2)+pow((Coord1.Getdata(i,3)-crd3.Getdata(0,3)),2));
		distance5=sqrt(pow((Coord1.Getdata(i,1)-crd3.Getdata(1,1)),2)+pow((Coord1.Getdata(i,2)-crd3.Getdata(1,2)),2)+pow((Coord1.Getdata(i,3)-crd3.Getdata(1,3)),2));
		distance6=sqrt(pow((Coord1.Getdata(i,1)-crd3.Getdata(2,1)),2)+pow((Coord1.Getdata(i,2)-crd3.Getdata(2,2)),2)+pow((Coord1.Getdata(i,3)-crd3.Getdata(2,3)),2));
		if ((fabs(distance4-edge4)<dem)&&(fabs(distance5-edge5)<dem)&&(fabs(distance6-edge6)<dem))
		{
			crd3.Setdata(3,Coord1.Getdata(i),Coord1.Getdata(i,1),Coord1.Getdata(i,2),Coord1.Getdata(i,3));

		}
	}
	for (int i=0;i<M;i++)
	{
		//AfxMessageBox(crd3.Getdata(i) ,MB_OK|MB_ICONINFORMATION);
	}
	

}




double CCordtrans::Calthreshold(CCrdNum crd1, CCrdNum crd2)//����ƥ����ֵ
{
	double Maxdist=0.0;
	double Dist=0.0;
	double threshold=3.0;
	int count1=crd1.GetmemNum();
	int count2=crd2.GetmemNum();
	if (count1==count2)
	{
		for (int i=0;i<count2;i++)
		{
			Dist=sqrt(pow((crd2.Getdata(i,1)-crd1.Getdata(i,1)),2)+pow((crd2.Getdata(i,2)-crd1.Getdata(i,2)),2)+pow((crd2.Getdata(i,2)-crd1.Getdata(i,2)),2));
			if (Dist>Maxdist)
			{
				Maxdist=Dist;
			}
		}
		threshold=5*Maxdist;
	} 
	else
	{
		AfxMessageBox("ͬ�����������",MB_OK|MB_ICONINFORMATION);
		
	}
	return threshold;
	
}

void CCordtrans::Finddistancepoint(int FIGURE,CCrdNum &crd3, CCrdNum &crd4)//Ѱ�Ҿ�������һ����ֵ�ĵ��
{
	int scount=0;
	CCrdNum crd30;
	CCrdNum crd31;
	CCrdNum crd40;
	CCrdNum crd1to2;
	CMatrix R0(3,3);
	CMatrix T0(3,1);
	double bound=3.0;
	//double RMS=0.0;
	double dist=0.0;
	double Mindist=0.0;
	CString tes;
	switch (FIGURE)
	{
	case 1:
		Findsamename(crd30,crd40);
		break;
	case 2:
		Distillpoint(crd30,crd40);
		break;
	}
	
	PrecisionCalPar(crd30,crd40,R0,T0);
 	Translate(R0,T0,crd30,crd31);
 	bound=Calthreshold(crd31,crd40);
	Translate(R0,T0,Coord1,crd1to2);

	for (int i=0;i<Ncount2;i++)
	{
		Mindist=sqrt(pow((Coord2.Getdata(i,1)-crd1to2.Getdata(0,1)),2)+pow((Coord2.Getdata(i,2)-crd1to2.Getdata(0,2)),2)+pow((Coord2.Getdata(i,3)-crd1to2.Getdata(0,3)),2));
		for (int j=1;j<Ncount1;j++)
		{
			dist=sqrt(pow((Coord2.Getdata(i,1)-crd1to2.Getdata(j,1)),2)+pow((Coord2.Getdata(i,2)-crd1to2.Getdata(j,2)),2)+pow((Coord2.Getdata(i,3)-crd1to2.Getdata(j,3)),2));
			if (dist<=Mindist)
			{
				Mindist=dist;
			}

		}
		if (Mindist<bound)
		{
			scount=scount+1;
			
		}
	}
	tes.Format("%d",scount);
	tes=tes+"����ƥ���ϵĵ���";
	AfxMessageBox(tes,MB_OK|MB_ICONINFORMATION);
	crd3.Init(scount);
	crd4.Init(scount);
	int k=0;
	int temp=0;
	for (int i=0;i<Ncount2;i++)
	{
		Mindist=sqrt(pow((Coord2.Getdata(i,1)-crd1to2.Getdata(0,1)),2)+pow((Coord2.Getdata(i,2)-crd1to2.Getdata(0,2)),2)+pow((Coord2.Getdata(i,3)-crd1to2.Getdata(0,3)),2));
		for (int j=0;j<Ncount1;j++)
		{
			dist=sqrt(pow((Coord2.Getdata(i,1)-crd1to2.Getdata(j,1)),2)+pow((Coord2.Getdata(i,2)-crd1to2.Getdata(j,2)),2)+pow((Coord2.Getdata(i,3)-crd1to2.Getdata(j,3)),2));
			if (dist<=Mindist)
			{
				Mindist=dist;
				temp=j;

			}

		}
		if (Mindist<bound)
		{
			crd3.Setdata(k,Coord1.Getdata(temp),Coord1.Getdata(temp,1),Coord1.Getdata(temp,2),Coord1.Getdata(temp,3));
			crd4.Setdata(k,Coord2.Getdata(i),Coord2.Getdata(i,1),Coord2.Getdata(i,2),Coord2.Getdata(i,3));
			k=k+1;	
		}
	}
	
}


CString  CCordtrans::Distancematch(CCrdNum &Coord1to2)//����ƥ��
{
	CString txt0="\n";
	CCrdNum crd3;
	CCrdNum crd4;
	CCrdNum crd1to2;
	CMatrix R0(3,3);
	CMatrix T0(3,1);
	Finddistancepoint(1,crd3,crd4);
	PrecisionCalPar(crd3,crd4,R0,T0);
	Translate(R0,T0,Coord1,crd1to2);
	Coord1to2=crd1to2;
	return Print(R0,T0,crd3,crd4,crd1to2);
}

CString CCordtrans::Distancematch(CMatrix &R, CMatrix &T, CCrdNum &Coord1to2)
{
	CString txt0="\n";
	CCrdNum crd3;
	CCrdNum crd4;
	CCrdNum crd1to2;
	CMatrix R0(3,3);
	CMatrix T0(3,1);
	Finddistancepoint(1,crd3,crd4);
	PrecisionCalPar(crd3,crd4,R0,T0);
	R=R0;
	T=T0;
	Translate(R0,T0,Coord1,crd1to2);
	Coord1to2=crd1to2;
	return Print(R0,T0,crd3,crd4,crd1to2);
	
}

CString  CCordtrans::Pointnamematch(CCrdNum &Coord1to2)//����ƥ�����
{
	CString txt0="\n";
	CCrdNum crd3;
	CCrdNum crd4;
	CCrdNum crd1to2;
	CMatrix R0(3,3);
	CMatrix T0(3,1);
	Findsamename(crd3, crd4);
	PrecisionCalPar(crd3,crd4,R0,T0);
	Translate(R0,T0,Coord1,crd1to2);
	Coord1to2=crd1to2;
	return Print(R0,T0,crd3,crd4,crd1to2);
}

CString CCordtrans::Pointnamematch(CMatrix &R, CMatrix &T, CCrdNum &Coord1to2)
{
	CString txt0="\n";
	CCrdNum crd3;
	CCrdNum crd4;
	CCrdNum crd1to2;
	CMatrix R0(3,3);
	CMatrix T0(3,1);
	Findsamename(crd3, crd4);
	PrecisionCalPar(crd3,crd4,R0,T0);
	//CalPar(crd3,crd4,R0,T0);
	R=R0;
	T=T0;
	Translate(R0,T0,Coord1,crd1to2);
	Coord1to2=crd1to2;
	return Print(R0,T0,crd3,crd4,crd1to2);
	
}


CString  CCordtrans::Initialmatch(CCrdNum &Coord1to2)//˳��ƥ��
{
	CString txt0="\n";
	CMatrix R0(3,3);
	CMatrix T0(3,1);
	CCrdNum crd1to2;
	PrecisionCalPar(R0,T0);
	Translate(R0,T0,Coord1,crd1to2);
	Coord1to2=crd1to2;
	return Print(R0,T0,Coord1,Coord2,crd1to2);
	
}

CString CCordtrans::Initialmatch(CMatrix &R, CMatrix &T, CCrdNum &Coord1to2)
{
	CString txt0="\n";
	CMatrix R0(3,3);
	CMatrix T0(3,1);
	CCrdNum crd1to2;
	PrecisionCalPar(R0,T0);
	//CalPar(R0,T0);
	R=R0;
	T=T0;
	Translate(R0,T0,Coord1,crd1to2);
	Coord1to2=crd1to2;
	return Print(R0,T0,Coord1,Coord2,crd1to2);
	
}


CString  CCordtrans::Freematch(CCrdNum &Coord1to2)//����ƥ��
{
	CString txt0="\n";
	CCrdNum crd3;
	CCrdNum crd4;
	CCrdNum crd1to2;
	CMatrix R0(3,3);
	CMatrix T0(3,1);
	Finddistancepoint(2,crd3,crd4);
	PrecisionCalPar(crd3,crd4,R0,T0);
	Translate(R0,T0,Coord1,crd1to2);
	Coord1to2=crd1to2;
	return Print(R0,T0,crd3,crd4,crd1to2);
}


CString CCordtrans::Freematch(CMatrix &R, CMatrix &T, CCrdNum &Coord1to2)
{
	CString txt0="\n";
	CCrdNum crd3;
	CCrdNum crd4;
	CCrdNum crd1to2;
	CMatrix R0(3,3);
	CMatrix T0(3,1);
	Finddistancepoint(2,crd3,crd4);
	PrecisionCalPar(crd3,crd4,R0,T0);
	R=R0;
	T=T0;
	Translate(R0,T0,Coord1,crd1to2);
	Coord1to2=crd1to2;
	return Print(R0,T0,crd3,crd4,crd1to2);

	
}



CString CCordtrans::Differencebetween(CCrdNum crd1, CCrdNum crd2)//�����λƫ��
{
	CString txt0="\n";
	int memsize1=crd1.GetmemNum();
	int memsize2=crd2.GetmemNum();
	CString detX,detY,detZ;
	double crdRMS=-1;
	double RMS=0;
	CString dets;
	CString CRMS;
	if (memsize1=memsize2)
	{
		txt0=txt0+"�ο���"+"  "+"������"+"     "+"��X"+"       "+"��Y"+"       "+"��Z"+"       "+"��s"+"\n";
		for (int i=0;i<memsize2;i++)
		{
			detX.Format("%0.4lf",crd2.Getdata(i,1)-crd1.Getdata(i,1));
			detY.Format("%0.4lf",crd2.Getdata(i,2)-crd1.Getdata(i,2));
			detZ.Format("%0.4lf",crd2.Getdata(i,3)-crd1.Getdata(i,3));
			crdRMS=sqrt(pow((crd2.Getdata(i,1)-crd1.Getdata(i,1)),2)+pow((crd2.Getdata(i,2)-crd1.Getdata(i,2)),2)+pow((crd2.Getdata(i,3)-crd1.Getdata(i,3)),2));
			dets.Format("%0.4lf",crdRMS);
			RMS=RMS+crdRMS*crdRMS;
			txt0=txt0+crd2.Getdata(i)+"      "+crd1.Getdata(i)+"      "+detX+"    "+detY+"    "+detZ+"    "+dets+"\n";
		}
		RMS=sqrt(RMS/memsize2);
		CRMS.Format("%0.4lf",RMS);
		txt0=txt0+"��λƫ��ľ�������RMS����";
		txt0=txt0+CRMS;

	}
	else
	{
		AfxMessageBox("����������㼯�е�ĸ��������",MB_OK|MB_ICONINFORMATION);
		txt0=txt0+"��������"+"\n";

	}
	return  txt0;


}

void CCordtrans::save()//�����ļ�
{
// 	CFileDialog dlg(FALSE,"txt"); // ����һ���ļ�����Ի������
// 	if(dlg.DoModal()==IDOK) 
// 	{
// 		CStdioFile File;
// 		File.Open(dlg.GetFileName(),CFile.modeCreate|CFile::modeWrite|CFile::modeNoTruncate);
// 		File.WriteString(txt);
// 		File.Close();
// 	}

}


//������ת��Ϊ��
//double radianΪ����Ļ���
double CCordtrans::r2d( double radian)
{
	double du,fen,miao;
	if (radian>0.0)
	{
		du=floor(radian*(180.0/PI));
		fen=floor((radian*(180.0/PI)-du)*60);
		miao=radian*(180.0/PI)*3600-du*3600-fen*60;
		return(du+fen/100.0+miao/10000.0);
	} 
	else
	{
		du=ceil(radian*(180.0/PI));
		fen=ceil((radian*(180.0/PI)-du)*60);
		miao=radian*(180.0/PI)*3600-du*3600-fen*60;
		return(du+fen/100.0+miao/10000.0);
	}
	
}


//���㷽���
//CMatrix R0Ϊ��ת����
CMatrix CCordtrans::DirectionAngle(CMatrix R0)
{
	CMatrix Angle(3,3);
	double radion=0.0;
	int NumRows=R0.GetNumRows();//��ȡ���������
	int NumColumns=R0.GetNumColumns();//��ȡ���������
	for (int i=0;i<NumRows;i++)
	{
		for(int j=0;j<NumColumns;j++)
		{
			radion=acos(R0.GetElement(i,j));
			Angle.SetElement(i,j,r2d(radion));
		}
	}
	return Angle;

}


CString CCordtrans::Print(CMatrix R0, CMatrix T0, CCrdNum crd3, CCrdNum crd4, CCrdNum crd1to2)
{
	CString txt0="\n";
	CCrdNum crd3to4;
	Translate(R0,T0,crd3,crd3to4);
	txt0=txt0+"ƽ�Ʋ�����"+"\n";
	txt0=txt0+T0.ToString("   ")+"\n";
	txt0=txt0+"\n";
	txt0=txt0+"��ת�ǣ�"+"\n";
	txt0=txt0+AngleofSwing(R0)+"\n";
	txt0=txt0+"\n";
	txt0=txt0+"��ת����"+"\n";
	txt0=txt0+R0.ToString("   ")+"\n";
	txt0=txt0+"\n";
	txt0=txt0+"����ת���ķ���ǣ�DEG����"+"\n";
	txt0=txt0+DirectionAngle(R0).ToString("   ");
	txt0=txt0+"\n"+"\n";
	txt0=txt0+"������ת�����ο�����ϵ�е����꣺"+"\n";
	txt0=txt0+crd1to2.Printpoint()+"\n";
	txt0=txt0+"��λƫ�RMS����"+"\n";
	txt0=txt0+Differencebetween(crd3to4,crd4);
	return txt0;	
}

CString CCordtrans::AngleofSwing(CMatrix R0, double Rx, double Ry, double Rz)//������ת��
{
	CString txt0="\n";
	CString crx,cry,crz;
	Rx=r2d(atan(-(R0.GetElement(1,2)/R0.GetElement(2,2))));
	Ry=r2d(asin(R0.GetElement(0,2)));
	Rz=r2d(atan(-(R0.GetElement(0,1)/R0.GetElement(0,0))));
	crx.Format("%0.8lf",Rx);//���������ת�Ƕ�
	cry.Format("%0.8lf",Ry);
	crz.Format("%0.8lf",Rz);
	txt0=txt0+"Rx="+crx+"\n"+"Ry="+cry+"\n"+"Rz="+crz+"\n";
	return txt0;
	
}

CString CCordtrans::AngleofSwing(CMatrix R0)//������ת��
{
	CString txt0="\n";
	CString crx,cry,crz;
	double Rx,Ry,Rz;
	Rx=r2d(atan(-(R0.GetElement(1,2)/R0.GetElement(2,2))));
	Ry=r2d(asin(R0.GetElement(0,2)));
	Rz=r2d(atan(-(R0.GetElement(0,1)/R0.GetElement(0,0))));
	crx.Format("%0.8lf",Rx);//���������ת�Ƕ�
	cry.Format("%0.8lf",Ry);
	crz.Format("%0.8lf",Rz);
	txt0=txt0+"Rx="+crx+"\n"+"Ry="+cry+"\n"+"Rz="+crz+"\n";
	return txt0;
	
}

void CCordtrans::AngleofSwing(CMatrix R0, CString &Rx, CString &Ry, CString  &Rz)//������ת��
{
	
	Rx.Format("%0.8lf",r2d(atan(-(R0.GetElement(1,2)/R0.GetElement(2,2)))));//���������ת�Ƕ�
	Ry.Format("%0.8lf",r2d(asin(R0.GetElement(0,2))));
	Rz.Format("%0.8lf",r2d(atan(-(R0.GetElement(0,1)/R0.GetElement(0,0)))));	
}

