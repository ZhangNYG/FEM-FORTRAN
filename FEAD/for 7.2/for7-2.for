C	MASTER TEMPERATURE
	REAL K0,K,NF
	INTEGER E,E1,E2,A,B,S,A0,A1,A2,TL,AB,E0
	LOGICAL AAA
	DIMENSION K(2000),N(0:300),K0(300),X(200),R(300),P(300)
     1,I(500),J(500),M(600),AL(200),TI(200),TL(30),T00(30)
	1,MB(0:20),AL0(20),AT0(20),KEY(2),MA(20),DQ(500),DQ0(500),G(500)
	COMMON X,R,P,K,C
c	************************************************
c	������Ʊ���key
	OPEN(1,FILE='dat.txt')
	OPEN(4,FILE='DAT1.txt',STATUS='NEW')
	READ(1,5)KEY
C	���key��1��=0 �������Գƺ㶨�¶ȳ����������ƽ���ȶ��¶ȳ�
	
	
	IF(KEY(1).EQ.0) WRITE(4,106)KEY(1),KEY(2)
	IF(KEY(1).NE.0) WRITE(4,110)KEY(1),KEY(2)
	
C	�����ˢ����ʼ����
	
	READ(1,8)L0,L01,E1,E2,A0,B,H,C
	READ(1,*)(X(L),L=1,L0)
	READ(1,*)(R(L),L=1,L0)
	READ(1,*)(I(E),E=1,E2)
	READ(1,*)(J(E),E=1,E2)
	READ(1,*)(M(E),E=1,E2)
	
	
	WRITE(4,1112)L0,L01,E1,E2,A0,B,H,C
	WRITE(4,1011)
1011	FORMAT(///,8X,'�ڵ��',3X,'X����',3X,'Y����')

	WRITE(4,1113)(L,X(L),R(L),L=1,L0)
	
	WRITE(4,1012)
1012	FORMAT(///,14X,'��Ԫ��',10X,'I',14X,'J',14X,'M',/)

	WRITE(4,2200)(E,I(E),J(E),M(E),E=1,E2)

	
C	������ڵ�һ��߽������������������¶�ֵΪ��֪�ı߽�ڵ�ź���Ӧ���¶�ֵ
	IF(L01.EQ.0)GOTO 101
               
	READ(1,11)TL,T00
	WRITE(4,1114)TL,T00
C	������ڵڶ���߽���������������������߽�Ԫ�����һ����Ԫ�ţ����鵥Ԫ�ڽ������¶Ⱥ��Ƚ���ϵ��

101	IF(B.EQ.0)GOTO 15
c*************************************************************�иĶ����鳬������
c************************************************************
	READ(1,12)(MB(III),III=1,20),AT0,AL0
c**************************************************************
c***********************************************************
	WRITE(4,1115)MB,AT0,AL0
	
C	������ڵڶ���߽��������򽫸���߽���Ƚ���ϵ������Χ�����¶ȷ��������ĸ�����Ԫ
	MB(0)=E1
13	DO 14 A=1,B
	DO 14 E=MB(A-1)+1,MB(A)
	L=E-E1
	TI(L)=AT0(A)
14	AL(L)=AL0(A)
15	CONTINUE
C	���ڵ����껹ԭ��ʵ������
	DO 2090 L=1,L0
	X(L)=X(L)/C
2090	R(L)=R(L)/C
C	�����ȴ���ϵ�����Ƚ���ϵ���ĵ�λ����
	H=H/1E2
	DO 220 A=1,E2-E1
220	AL(A)=AL(A)/1E4
C	�����¶ȸնȾ���k���ĶԽ���Ԫ����һά������������еĵ�ַ�����һ���Խ���Ԫ�صĵ�ֵַ
C	�����¶ȸնȾ���һά����洢ʱ��Ԫ������
	N(0)=0
	DO 24 L=1,L0
	A=0
	DO 22 E=1,E2
	AAA=((I(E).EQ.L).OR.(J(E).EQ.L).OR.(M(E).EQ.L))
c	�����д˵�Ԫ��һ�ڵ����Lʱ���ɽ�����һ��
	IF(.NOT.AAA)GOTO 22
	A1=I(E)
	A2=L-A1
	IF(A2.GT.A)A=A2
	A1=J(E)
	A2=L-A1
	IF(A2.GT.A)A=A2
	A1=M(E)
	A2=L-A1
	IF(A2.GT.A)A=A2
22	CONTINUE
	N(L)=N(L-1)+A+1
24	CONTINUE
	S=N(L0)
	WRITE(4,30)S
C	���¶ȸնȾ��󷽳����Ҷ������
40	DO 50 A=1,S
50	K(A)=0.0
	DO 60 L=1,L0
60	P(L)=0.0
C	����ϵ��bi��bj��bm��ci��cj��cm�͡�e��
	DO 70 E=1,E2
	IE=I(E)
	JE=J(E)
	ME=M(E)
	XI=X(IE)
	XJ=X(JE)
	XM=X(ME)
	RI=R(IE)
	RJ=R(JE)
	RM=R(ME)
	BI=RJ-RM
	BJ=RM-RI
	BM=RI-RJ
	CI=XM-XJ
	CJ=XI-XM
	CM=XJ-XI
	DELTA=(BI*CJ-BJ*CI)/2
	SI=SQRT(BI*BI+CI*CI)
C	�������ݴ��������Ԫ���С�ڻ�����㣬���ӡ�õ�Ԫ��ţ���ͣ��
	IF(DELTA)201,201,202
201	WRITE(4,222)E
	STOP 222
202	CONTINUE
c	����Ǳ߽絥Ԫ���� ����
	IF(E.GT.E1)DI=AL(E-E1)
	IF(E.GT.E1)FI=DI*TI(E-E1)
C	�������ƽ���ȶ��¶ȳ�����ri=rj=rm=1
	IF(KEY(1).EQ.0)GOTO 65
	RI=1.0
	RJ=1.0
	RM=1.0
65	AA=H*(RI+RJ+RM)/12./DELTA
C
C*********************************
C	���㵥Ԫ�¶ȸնȾ����������ξ�������Ԫ�أ���������ӵ������¶ȸնȾ�����
	A1=N(IE)
	K(A1)=K(A1)+AA*SI*SI
	A1=N(JE)
	IF(E.LE.E1)
     1K(A1)=K(A1)+AA*(BJ*BJ+CJ*CJ)
	IF(E.GT.E1)
     1K(A1)=K(A1)+AA*(BJ*BJ+CJ*CJ)+DI/4.*SI*(RJ+RM/3.)
	A1=N(ME)
      IF(E.LE.E1)
	1K(A1)=K(A1)+AA*(BM*BM+CM*CM)
	IF(E.GT.E1)
     1K(A1)=K(A1)+AA*(BM*BM+CM*CM)+DI/4.*SI*(RM+RJ/3.)
	IF(IE.GT.JE)A1=N(IE)-(IE-JE)
	IF(IE.LE.JE)A1=N(JE)-(JE-IE)
	K(A1)=K(A1)+AA*(BI*BJ+CI*CJ)
	IF(JE.GT.ME)A1=N(JE)-(JE-ME)
	IF(JE.LE.ME)A1=N(ME)-(ME-JE)
	IF(E.LE.E1)K(A1)=K(A1)+AA*(BJ*BM+CJ*CM)
	IF(E.GT.E1)
     1K(A1)=K(A1)+AA*(BJ*BM+CJ*CM)+DI/12.*SI*(RJ+RM)
	IF(IE.GT.ME)A1=N(IE)-(IE-ME)
	IF(IE.LE.ME)A1=N(ME)-(ME-IE)
	K(A1)=K(A1)+AA*(BI*BM+CI*CM)
C	�����Ҷ���
	IF(E.GT.E1)P(IE)=P(JE)+FI/3.*SI*(RJ+RM/2.)
	IF(E.GT.E1)P(ME)=P(ME)+FI/3.*SI*(RM+RJ/2.)
70	CONTINUE
C	�����֪�ڵ��¶ȱ߽�����������նȾ��󣬼����Խ���Ԫ�س�һ�������������¶�ֵ������Ӧ���Ҷ���
	IF(L01.EQ.0)GOTO 90
	DO 205 A=1,L01
	L=TL(A)
	AA=T00(A)
	A1=N(L)
	K(A1)=K(A1)*1E8
	P(L)=K(A1)*AA
205	CONTINUE
90	CONTINUE
C		������Դ��������飬����ڵ��¶�ֵ
	CALL PG2Y(L0,S,K,N,P)
	WRITE(4,121)(L,P(L),L=1,L0)
C	���Ҫ��������ߣ�������Ӧ��������������������
	IF(KEY(2).EQ.0)GOTO 990
	READ(1,131)T0,T1,DT,E0
	IF(A0.GT.0)READ(1,132)MA
C	����߽��ϵ����������ӵ�һ���߽絥Ԫ��ʼѭ�����ȼ���ϵ��bi��ci��si��
	CLOSE(1)
130	NF=0
	IEE=E0+1
	IE=I(E)
	JE=J(E)
	ME=M(E)
	XI=X(IE)
	XJ=X(JE)
	XM=X(ME)
	RI=R(IE)
	RJ=R(JE)
	RM=R(ME)
	BI=RJ-RM
	CI=XM-XJ
	SI=SQRT(BI*BI+CI*CI)
C	�������ƽ���ȶ��¶ȳ�������AA=1�����������Գ��ȶ��¶ȳ�������AA=��R1+Rm����
c	��EΪ��һ��߽絥Ԫʱ������Bj��Bm��Cj��Cm���͡�e
	A=E-E0
	IF(KEY(1).EQ.0)AA=(RJ+RM)*3.1415926
	IF(E.GT.E1)GOTO 160
	BJ=RM-RI
	BM=RI-RJ
	CJ=XI-XM
	CM=XJ-XI
	DELTA=(BI*CJ-BJ*CI)/2
C	����߽絥Ԫ�߽�jm�ϵĵ�Ԫʱ����������������������������е���
	DQ(A)=-H/2./DELTA*(SI*SI*P(IE)+(BI*BJ+CI*CJ)*P(JE)+(BI*BM
     1+CI*CM)*P(ME))*AA
	NF=NF+DQ(A)
	GOTO 170
C	����ǵڶ���߽絥Ԫ�����㵥λʱ��������������߽�jm���������������е���
160	A1=E-E1
	FI=TI(A1)
	DI=AL(A1)
	DQ(A)=DI*(FI-(P(JE)+P(ME))/2)*AA*SI
	NF=NF+DQ(A)
170	CONTINUE
C	����߽��ϸ���ÿСʱ�����������������
	MA(0)=E0
	AB=A0+B
	DO 200 A=1,AB
	DQ0(A)=0.0
C	���ڵ�һ��߽絥Ԫ������A0���ϵ�������
	IF(A.GT.A0)GOTO 190
	MMA=MA(A-1)+1
	MAA=MA(A)
	DO 180 E=MMA,MAA
	L=E-E0
	DQ(A)=DQ0(A)+DQ(L)
180	CONTINUE
C	���ڵڶ���߽絥Ԫ������B���ϵ�������
190	MMB=MB(A-A0-1)+1
	MBB=MB(A-A0)
	DO 200 E=MMB,MBB
	L=E-E0
	DQ0(A)=DQ0(A)+DQ(L)
200	CONTINUE
C	����߽����������������߽�ͱ߽絥Ԫjm���ϵ�������
	WRITE(4,300)NF,(LI,DQ0(LI),LI=1,AB),(LI,DQ(LI),LI=1,E2-E0)
C	����������ϸ�����������Ӧ��Ԫ��
	AA=T0
350	A=1
	DO 700 E=1,E2
	IE=I(E)
	JE=J(E)
	ME=M(E)
	IF((P(IE)-AA)*(P(JE)-AA).GT.0)GOTO 500	
400	I0=JE
	J0=ME
	CALL TT(I0,J0,A,AA,E)
500	IF((P(JE)-AA)*(P(ME)-AA).GT.0)GOTO 600
	I0=ME
	J0=ME
	CALL TT(I0,J0,A,AA,E)
600	IF((P(ME)-AA)*(P(IE)-AA).GT.0)GOTO 700
	I0=ME
	J0=IE
	CALL TT(I0,J0,A,AA,E)
700	CONTINUE
800	K(A)=0
	IF(A.LE.1)GOTO 900
	DO 810 L=1,A
810	G(L)=K(L)
C	��������ߵ��¶�ֵ����Ԫ�ź͸�������
	WRITE(4,820)AA,(G(L),L=1,A)
900	CONTINUE
C	���δ����������µ����ߣ�������������
	AA=AA+DT
	IF(AA.LE.T1)GOTO 350
				
C	*******************************************
C	��ʽ���
5	FORMAT(2I5)
106	FORMAT(15X,'* AXIS TEMPERATURE FIELD *'///,I5,I5)
110	FORMAT(15X,'* PLANE TEMPERATURE FILED *'///,I5,I5)
8	FORMAT(6I4,2F10.4)
1111	FORMAT(10F8.2)
2222	FORMAT(10I7)
1112	FORMAT(1X,'L0=',I2,5X,'L01=',I2,5X,'E1=',I2,5X,'E2=
     1',I2,5X,'A0=',I2,5X,'B=',I2,5X,'H=',F10.3,5X,'C=',F5.3)
1113	FORMAT(I10,2F10.3)
2200	FORMAT(4I15)
1114	FORMAT(////'TL:',30I4//,'T00:',30F10.2///)
11	FORMAT(30I4/30F10.2)
1115	FORMAT(////'MB:',(21I4)///,'AT0:',(20F10.4)///,'AL0:',(20F10.4)///)
12	FORMAT(20I4/20F8.1/20F10.4)
30	FORMAT(///10X,'2HS=',I5///)
222	FORMAT(10X,2HE,I5)
121	FORMAT(/10X,'TEMPERATURE OF NODAL POINTS'//(8F8.2))
131	FORMAT(3F10.3,I10)
132	FORMAT(20I4)
300	FORMAT(///10X,'NF=',F10.2///(I,F8.2)///(I,F8.2))
820	FORMAT(///10X,'T=',F10.2/3X,'E',10X,'X',10X,'R',/(3F10.3))
990	CONTINUE
	CLOSE(4)
	STOP
	END
C	��������ϵĵ���������Ӧ��Ԫ�ŵ��ӳ���
	SUBROUTINE TT(I0,J0,A,AA,E)
	INTEGER A,E
	REAL K
	DIMENSION X(300),R(300),P(300),K(2000)
	COMMON X,R,P,K,C
C	����������ϵĵ��x����
	XI=X(I0)
	RI=R(I0)
	XJ=X(J0)
	RJ=R(J0)
	AAL=(P(I0)-AA)/(P(I0)-P(J0))
	K(A)=E
	A=A+1
	K(A)=((1-AAL)*XI+AAL*XJ)*C
C	����������ϵ��r����
	A=A+1
	K(A)=((1-AAL)*RI+AAL*RJ)*C
	A=A+1
	K(A)=0
	A=A+1
	RETURN
	END
C*************************************************
C	�ӳ���������Դ����������ӳ���
	SUBROUTINE PG2Y(N,NQ,SK,NA,RF)	
	DIMENSION NA(N),RF(N),SK(NQ)
	DO 1 I=1,N
	II=I-1
	ID=NA(I)
	IF(I.LE.1)GOTO 1
	MI=I-NA(I)+NA(II)+1
	DO 2 J=MI,I
	IGP=NA(I)+J-I
	JD=NA(J)
	J1=J-1
	IF(J1.LE.0)GOTO 3
	MJ=J-NA(J)+NA(J1)+1
	IJ=MJ
	IF(MI.GT.MJ)IJ=MI
	IF(IJ.GT.J1)GOTO 3
	DO 4 K=IJ,J1
	K1=NA(I)+K-I
	K2=NA(J)+K-J
	KD=NA(K)
4	SK(IGP)=SK(IGP)-SK(K)*SK(K2)*SK(KD)
3	IF(J.EQ.I)GOTO 1
	SK(IGP)=SK(IGP)/SK(JD)
2	RF(I)=RF(I)-SK(IGP)*SK(JD)*RF(J)
1	RF(I)=RF(I)/SK(ID)
	NX=N-1
	DO 5 I=1,NX
	NXR=N-I
	NXD=NXR+1
	DO 5 J=NXD,N
	LR=J+1-NA(J)+NA(J-1)
	IF(LR.GT.NXR)GOTO 5
	LK=NA(J)-J+NXR
	RF(NXR)=RF(NXR)-SK(LK)*RF(J)
5	CONTINUE
	RETURN
	END											

					
			
				


	
		














	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		     