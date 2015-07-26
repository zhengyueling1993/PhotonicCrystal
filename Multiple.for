c     这个程序可以计算有限二维光子晶体的透射谱。
        IMPLICIT NONE    
c	 数据说明：
c      k0波矢，radium圆柱的半径，lattice周期长度。
c      Hank第三类Bessel函数 Hankr第三类Bessel函数的导数                                                           
c      JN第一类Bessel函数  JNr第一类Bessel函数的导数
c	 bessi修正的第一类Bessel函数,rbessi修正的第一类Bessel函数的导数
c	 x(N),Y(N)代表构成光子晶体的每个柱子的位置的X,Y坐标。N是总的柱子数
c	 D代表边界条件导致的系数，详细的推导可在我的论文中找到
c	 B是每个柱子的不同阶散射系数满足的线性方程组的左边。
c	 A是入射外场在每个柱子的不同阶展开系数。
C      XXX是每个柱子的不同阶散射系数
C      freq是频率，它和波矢K0满足如下关系freq=c*k0
c	 k是柱子里的波矢它和k0满足如下关系k=sqrt(epsi)*k0
c      epsi是柱子的介电常数。
c	 matrr存放散射系数满足的线性方程组左边的B的实部，
C	 matri存放散射系数满足的线性方程组左边的B的虚部，
C	 BR存放散射系数满足的线性方程组右边的A的实部，解完线性方程组以后存放的是
c	 散射系数的实部
c	 Bi存放散射系数满足的线性方程组右边的A的实部，解完线性方程组以后存放的是
c	 散射系数的虚部

        REAL*8 K0,RADIUM,LATTICE,PI                                 
        REAL*8 X1,Y1,X2,Y2,XX,YY
        COMPLEX*16 HANK,HANKR                                                          
        real*8 JN,JNR,bessi,rbessi, x(145),y(145)                                            
        REAL*8 ANGLE,R,temp1,k                                        
        COMPLEX*16 D(-2:2),B(725,725),A(725),xxx(725)                          
        COMPLEX*16 temp2,temp3,temp4,TEMP8                                      
        real*8  freq, RESULT,result1                                           
        REAL*8 MATRR(725,725),MATRI(725,725),BR(725,1),BI(725,1)
        REAL*8 JS(725),anglee                                     
        INTEGER INT1,INT2,INT3,INT4,M,M1,M2,L
        real*8 epsi
c     文件2调用的是每个柱子的XY平面上的坐标   
        OPEN (UNIT=2,FILE="e:\chenbin\pbg\circle1.dat")
        OPEN (UNIT=3,FILE="e:\chenbin\pbg\test1.dat")

        pi=3.1415926                                                                   
        LATTICE=1.0                                                                   
        RADIUM=0.45                                  
        EPSI=5.0
        
        do 1 int1=1,145
        read(2,*) x(int1),y(int1)
          x1=x(int1)
          x(int1)=x1+lattice/2.0                                               
1       CONTINUE                                                                      


        DO 2 freq=0.03, 0.7, 0.002
c	  频率是以2*Pi*C/a做为约化频率
	  PRINT*,"Please don't close this program"                                                                          
c	 k0：背景的波矢量， K:柱子的波矢量 
        k=freq*2.0*pi/lattice
        k0=k*sqrt(epsi) 

c      计算D矩阵的矩阵元，如果每个柱子的半径和介电常数一样那么这将对与
c      每个柱子都是一样的。                                                             
        DO 3 M=-2,2
        IF (M.LT.0) THEN                                                               
        TEMP1=K0*JNR(ABS(M),K0*RADIUM)*JN(ABS(M),K*RADIUM)                              
     $         -K*JN(ABS(M),K0*RADIUM)*JNR(ABS(M),K*RADIUM)                       
        TEMP2=K*JNR(ABS(M),K*RADIUM)*HANK(ABS(M),K0*RADIUM)                           
     $        -K0*HANKR(ABS(M),K0*RADIUM)*JN(ABS(M),K*RADIUM)                      
        ELSE                                                                           
        TEMP1=K0*JNR(M,K0*RADIUM)*JN(M,K*RADIUM)                                        
     $        -K*JN(M,K0*RADIUM)*JNR(M,K*RADIUM)                                 
        TEMP2=K*JNR(M,K*RADIUM)*HANK(M,K0*RADIUM)                                     
     $        -K0*HANKR(M,K0*RADIUM)*JN(ABS(M),K*RADIUM)                                   
        END IF                                                                         
        D(M)=TEMP1/TEMP2                                                             
3       CONTINUE                                                                      


c     设置B矩阵。                                             
        INT3=0                                                                         
        DO 4 INT1=1,145
        X1=X(INT1)                                                                     
        Y1=Y(INT1)                                                                     
        DO 4 M1=-2,2                                                                   
        INT3=INT3+1                                                                    
        INT4=0                                                                         
                                                                                
        DO 4 INT2=1,145                                                                
        X2=X(INT2)                                                                     
        Y2=Y(INT2)                                                                     
        DO 5 M2=-2,2                                                                   
        INT4=INT4+1                                                                    
                                                                                
        IF (INT1.EQ.INT2) THEN
         B(INT3,INT4+2+M1)=-1.0/D(M1)                                                 
         INT4=INT4+4 
          goto 4
        END IF                                                                         
                                                                                
        M=M2-M1                                                                        
        XX=X2-X1                                                                       
        YY=Y2-Y1                                                                       
        CALL TRAN(XX,YY,ANGLE,r)                                                       
                                                                                
        IF (M.LT.0) THEN                                                               
          TEMP4=(HANK(ABS(M),K0*r))*((-1)**M )                                    
        ELSE                                                                           
          TEMP4=HANK(M,K0*r)                                                             
        END IF                                                                         
c      
c  	 angle 是第i个柱子和第j个柱子柱心连线矢量在极坐标下的角分量

          TEMP3=CMPLX(0.0,M*(ANGLE+pi) )                                                 
          B(INT3,INT4)=EXP(TEMP3)*TEMP4                                                  
        CONTINUE                                                                      
4       CONTINUE                                                                      

c      计算外场在每个柱子的不同阶的展开系数。这里的外场
c      是一个在原点发射出的Bessel函数形式的外场。                                                                     
        INT3=0
      DO 6 INT1=1,145                                                                
       X1=X(INT1)                                                                     
       Y1=Y(INT1)                                                                     
       CALL TRAN(x1,Y1,angle,r)  
        DO 6 M1=-2,2                                                             
           INT3=INT3+1
      IF (M1.LT.0) THEN
        TEMP8=HANK( ABS(M1),K0*R )*((-1)**M1)
        ELSE
        TEMP8=HANK(M1,K0*R)
        END IF
        TEMP2=CMPLX( 0.0,-M1*(ANGLE+PI) )
        A(INT3)=TEMP8*EXP(TEMP2)
6     CONTINUE


        DO 9 INT1=1,725
        BR(INT1,1)=DREAL(-a(int1))  
          BI(INT1,1)=DIMAG(-a(int1))                                               
          DO 9 INT2=1,725
          MATRR(INT1,INT2)=DREAL(B(INT1,INT2))                                                                     
        MATRI(INT1,INT2)=DIMAG(B(INT1,INT2))
9       CONTINUE

        L=0
        CALL ACJDN(MATRR,MATRI,725,BR,BI,1,L,JS)                                                                          
c     这是调用计算线性方程组的程序， 

        DO 44 INT1=1,725
        XXX(INT1)=CMPLX( BR(INT1,1),BI(INT1,1) )
44      CONTINUE

c      对不同角度的透射系数积分将得到全能带的信息.
对
          result1=0.0
        do 1111 anglee=0.0,2.0*pi,0.03 
        INT2=0  

        TEMP2=CMPLX(0.0, 0.0)
        TEMP4=CMPLX(0.0, -1.0)   
                                                                               
        DO 11 INT1=1,145
        X1=X(INT1)
        Y1=Y(INT1)
        CALL TRAN(X1,Y1,ANGLE,R)
        DO 11 M=-2,2
        TEMP3=CMPLX(0.0, -K0*R*COS(ANGLE-anglee)+m*anglee )
        INT2=INT2+1
        TEMP2=TEMP2+EXP(TEMP3)*(TEMP4**M)*XXX(INT2)
11      CONTINUE

        result=2.0*( (abs(1+temp2))**2 )/( pi*k0 )
      result1=result1+result*0.03
1111  continue

        WRITE(3,*) FREQ, log( RESULT1 )
        print*, freq,log( result1)

2       CONTINUE
        END                                                                     


c	 该子程序将直角平面的坐标改成极坐标的表示。 
        subroutine tran(x,y,angle,r)
        real*8 x,y,angle,r,x0,y0,temp,pi
        pi=3.1415926
        r=sqrt(x*x+y*y)
        x0=0.0
        y0=0.0
        if (x. GT.x0) then
           if (y.GT.y0) then
           temp=y/x
           angle=atan(temp)
           else if (y.EQ.y0) then
           angle=0.0
           else
           temp=abs(y/x)
           angle=2.0*pi-atan(temp) 
           end if
        else if (x.EQ.x0) then
           if (y.GT.y0) then
           angle=pi/2.0
           else if (y.EQ.y0) then
           angle=0.0
           else
           angle=3.0*pi/2.0
           end if
        else
           if (y.GT.y0) then
           temp=abs(y/x)
           angle=pi-atan(temp)
           else if (y.EQ.y0) then
            angle=pi
           else
            temp=abs(y/x)
                angle=pi+atan(temp)
           end if
        end if
        return
        end               
                                                                                             
c       该子程序计算第三类Bessel函数的导数                                                                         
        FUNCTION HANKR(N,X)                                                            
        INTEGER N                                                                      
        REAL*8 X                                                                       
        COMPLEX*16 HANKR,HANK                                                            
        IF(N.EQ.0) THEN                                                                
        HANKR=-HANK(1,X)                                                               
        ELSE                                                                           
        HANKR=-N*HANK(N,X)/X+HANK(N-1,X)                                               
        END IF                                                                         
        RETURN                                                                         
        END                                                                            

c     该子程序计算第三类Bessel函数                                                                                
        FUNCTION HANK(N,X)                                                             
        INTEGER N                                                                      
        REAL*8 X                                                                       
        COMPLEX*16 HANK                                                                
        REAL*8 JN,NN,J0,N0,J1,N1                                                       
        IF(N.GE.2) THEN                                                                
        HANK=CMPLX(JN(N,X),NN(N,X))                                                    
        ELSE IF (N.EQ.0) THEN                                                          
        HANK=CMPLX(J0(X),N0(X))                                                        
        ELSE                                                                           
        HANK=CMPLX(J1(X),N1(X))                                                        
        END IF                                                                         
        RETURN                                                                         
        END                                                                            
      
c	 该子程序计算第一类Bessel函数的导数                                                                                 
        FUNCTION JNR(N,X)                                                              
        INTEGER N                                                                      
        REAL*8 X,JNR,JN                                                                
        IF (N.EQ.0) THEN                                                               
        JNR=-JN(1,X)                                                                   
        ELSE                                                                           
        JNR=-N*JN(N,X)/X+JN(N-1,X)                                                     
        END IF                                                                         
        RETURN                                                                         
        END                                                                            
                                                                                
c	  该子程序计算第一类Bessel函数
        FUNCTION JN(N,X)                                                               
        INTEGER N,J                                                                    
        REAL*8 X,TOX,J1,J0,JN,BYP,BYM,BY                                               
                                                                                
        IF (X.EQ.0.0) THEN                                                             
        IF (N.EQ.0) THEN                                                               
        JN=1.0                                                                         
        ELSE                                                                           
        JN=0.0                                                                         
        END IF                                                                         
        GOTO 99                                                                        
        END IF                                                                         
                                                                                
        IF (N.GE.2) THEN                                                               
        TOX=2.0/X                                                                      
        BY=J1(X)                                                                       
        BYM=J0(X)                                                                      
        DO 11 J=1,N-1                                                                  
        BYP=J*TOX*BY-BYM                                                               
        BYM=BY                                                                         
        BY=BYP                                                                         
11      CONTINUE                                                                     
        JN=BY                                                                          
        ELSE IF (N.EQ.0) THEN                                                          
        JN=J0(X)                                                                       
        ELSE                                                                           
        JN=J1(X)                                                                       
        END IF                                                                         
99      RETURN                                                                       
        END                                                                            

c      该子程序计算第二类Bessel函数                                                                                
        FUNCTION NN(N,X)                                                               
        INTEGER N                                                                      
        REAL*8 X,TOX,N1,N0,NN,BYP,BYM,BY                                               
        IF (N.GE.2) THEN                                                               
        TOX=2.0/X                                                                      
        BY=N1(X)                                                                       
        BYM=N0(X)                                                                      
        DO 12 J=1,N-1                                                                  
        BYP=J*TOX*BY-BYM                                                               
        BYM=BY                                                                         
        BY=BYP                                                                         
12      CONTINUE                                                                     
        NN=BY                                                                          
        ELSE IF(N.EQ.0) THEN                                                           
        NN=N0(X)                                                                       
        ELSE                                                                           
        NN=N1(X)                                                                       
        END IF                                                                         
        RETURN                                                                         
        END                                                                            
                                                                                
c      该子程序计算第一类Bessel函数的零阶函数                                                                          
        FUNCTION J0(X)                                                                 
        REAL*8 X,J0                                                                    
        REAL*8 T,PI,P0,Q0,TEMP                                                         
        PI=3.1415926                                                                   
      IF (X.LE.4.0) THEN                                                        
      T=X/4.0                                                                   
      J0=1-3.9999998721*T**2+3.9999973021*T**4-1.7777560599*T**6                
     $ +0.4443584263*T**8-0.0709253492*T**10+0.0076771853*T**12                 
     $ -0.0005014415*T**14                                                      
        ELSE                                                                           
        T=4.0/X                                                                        
      P0=SQRT(2*PI)*(0.3989422793-0.0017530620*T**2+0.0001734300*T**4           
     $ -0.0000487613*T**6+0.0000173565*T**8-0.0000037043*T**10)                 
      Q0=T*sqrt(2*PI)*(-0.0124669441+0.0004564324*T**2                          
     $ -0.0000869791*T**4+0.0000342468*T**6-0.0000142078*T**8                   
     $ +0.0000032312*T**10)                                                     
      TEMP=SQRT(2.0/(PI*X))                                                     
      J0=TEMP*(P0*COS(X-PI/4.0)-Q0*SIN(X-PI/4.0))                               
        END IF                                                                         
        RETURN                                                                         
        END               
	                                                               
c      该子程序计算第一类Bessel函数的一阶函数                                                                                 
        FUNCTION J1(X)                                                                 
        REAL*8 X,J1                                                                    
        REAL*8 PI,T                                                                    
        PI=3.1415926                                                                   
      IF(X.LE.4.0) THEN                                                         
      T=X/4.0                                                                   
      J1=(1.9999999998-3.9999999710*T**2+2.6666660544*T**4                      
     $ -0.8888839649*T**6+0.1777582922*T**8-0.0236616773*T**10                  
     $ +0.0022069155*T**12-0.0001289769*T**14)*T                                
        ELSE                                                                           
      T=4.0/X                                                                   
      P1=sqrt(2*PI)*(0.3989422819+0.0029218256*T**2-0.0002232030*T**4           
     $ +0.0000580759*T**6-0.0000200920*T**8+0.0000042414*T**10)                 
      Q1=t*sqrt(2*PI)*(0.0374008364-0.0006390400*T**2+0.0001064741*T**4         
     $ -0.0000398708*T**6+0.0000162200*T**8-0.0000036594*T**10)                        
        TEMP=SQRT(2.0/(PI*X))                                                          
      J1=TEMP*(P1*COS(X-3.0*PI/4.0)-Q1*SIN(X-3.0*PI/4.0))                       
        END IF                                                                         
        RETURN                                                                         
        END                                                                            

c      该子程序计算二第类Bessel函数的零阶函数                                                                                 
        FUNCTION N0(X)                                                                 
        REAL*8 X,N0                                                                    
        REAL*8 TEMP1,PI,T,J0,P0,Q0,TEMP2                                               
        PI=3.1415926                                                                   
        IF (X.LE.4.0) THEN                                                             
        T=X/4.0                                                                        
      TEMP1=-0.5772156649-1.6911374142*T**2+3.6911388793*T**4                   
     $ -2.2331102234*T**6+0.6694321484*T**8-0.1214187561*T**10                  
     $ +0.0148999271*T**12-0.0013508487*T**14+0.0000891322*T**16                
      N0=2.0*(J0(X)*log(X/2.0)-TEMP1)/(PI)                                      
        ELSE                                                                           
      T=4.0/X                                                                   
      P0=sqrt(2*pi)*(0.3989422793-0.0017530620*t**2+0.0001734300*t**4           
     $ -0.0000487613*t**6+0.0000173565*t**8-0.0000037043*t**10)                 
      Q0=t*sqrt(2*pi)*(-0.0124669441+0.0004564324*t**2                          
     $ -0.0000869791*t**4+0.0000342468*t**6-0.0000142078*t**8                   
     $ +0.0000032312*t**10)                                                     
      TEMP2=SQRT(2.0/(PI*X))                                                    
      N0=TEMP2*(P0*SIN(X-PI/4.0)+Q0*COS(X-PI/4.0))                              
        END IF                                                                         
        RETURN                                                                         
        END                                                                            

c      该子程序计算第二类Bessel函数的一阶函数                                                                                 
        FUNCTION N1(X)                                                                 
        REAL*8 N1,X                                                                    
        REAL*8 TEMP2,T,PI,J1,TEMP3                                                     
        PI=3.1415926                                                                   
        IF (X.LE.4.0) THEN                                                             
        T=X/4.0                                                                        
      TEMP2=(1.0000000004-0.6177253972*T**2-10.7645472724*T**4                  
     $ +11.6207891416*T**6-4.9105291148*T**8+1.1418033012*T**10                 
     $ -0.1691081720*T**12+0.0169921876*T**14-0.0010266368*T**16)/(4*T)         
      N1=2.0*(J1(X)*log(X/2.0)-TEMP2)/PI                                        
        ELSE                                                                           
      T=4.0/X                                                                   
      P1=SQRT(2*PI)*(0.3989422819+0.0029218256*T**2-0.0002232030*T**4           
     $ +0.0000580759*T**6-0.0000200920*T**8+0.0000042414*T**10)                 
      Q1=T*sqrt(2*PI)*(0.0374008364-0.0006390400*T**2+0.0001064741*T**4         
     $ -0.0000398708*T**6+0.0000162200*T**8-0.0000036594*T**10)                        
      TEMP3=SQRT(2.0/(PI*X))                                                    
      N1=TEMP3*(P1*SIN(X-3.0*PI/4.0)+Q1*COS(X-3.0*PI/4.0))                      
        END IF                                                                         
        RETURN                                                                         
        END                                                                            

c      该子程序计算第一类修正Bessel函数的导数         
        function rbessi(n,x)
        real*8 bessi,x,rbessi
        integer n
        IF (N.EQ.0) THEN
        rbessi=-bessi(1,x)
        ELSE
        rbessi=-n*bessi(n,x)/x+bessi(n-1,x)
        END IF                                                                         
        RETURN
        end

c      该子程序计算第一类修正Bessel函数 
        function bessi(n,x)
        integer n
        real*8 x,bigno,bigni,bessi0,bessi1
        real*8 bim,bip,bi,bessi,tox
        integer iacc,m,j
        parameter (iacc=40,bigno=1.0e10,bigni=1.0e-10)
        if (n.LT.2) then
          if (n.EQ.0) then
          bessi=bessi0(x)
          else
          bessi=bessi1(x)
          end if
        else
        tox=2.0/x
        bip=0.0
        bi=1.0
        bessi=0.0
        m=2*((n+int(sqrt(float(iacc*n)))))
        do 111 j=m,1,-1
           bim=bip+float(j)*tox*bi
           bip=bi
           bi=bim
          if (abs(bi).GT.bigno) then
           bessi=bessi*bigni
           bi=bi*bigni
           bip=bip*bigni
           end if
         if (j.EQ.n) then
          bessi=bip
         end if
111     continue  
        bessi=bessi*bessi0(x)/bi
        end if
        return
        end

c      该子程序计算第一类修正Bessel函数的零阶函数 
        function bessi0(x)
        real*8  p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
        real*8  x,bessi0,y,ax
        data p1,p2,p3,p4,p5,p6,p7 /1.0d0,3.5156229d0,3.0899424d0,
     $    1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
        data q1,q2,q3,q4,q5,q6,q7,q8,q9 /0.39894228d0, 0.1328592d-1,
     $  0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,
     $  0.2635537d-1,-0.1647633d-1,0.392377d-2 /
        if (abs(x).LT.3.75) then
          y=(x/3.75)**2
          bessi0=p1+y*(p2+ y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
        else
          ax=abs(x)
          y=3.75/ax
          bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3
     $  +y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
        end if
        return
        end

c      该子程序计算第一类修正Bessel函数的一阶函数
        function bessi1(x)
        real*8   p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
        real*8   bessi1,y,x,ax
        data p1,p2,p3,p4,p5,p6,p7 /0.5d0, 0.87890594d0,0.51498869d0,
     $  0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3 /
        data q1,q2,q3,q4,q5,q6,q7,q8,q9 / 0.39894228d0,-0.3988024d-1,
     $-0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1,
     $0.1787654d-1,-0.420059d-2 /
        if (abs(x).LT.3.75) then
          y=(x/3.75)**2
          bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
        else
          ax=abs(x)
          y=3.75/ax
          bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4
     $  +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
        end if
        return
        end


c	  该子程序计算线性方程组。
      SUBROUTINE ACJDN(AR,AI,N,BR,BI,M,L,JS)
        INTEGER N,M,L
        REAL*8 AR(N,N),AI(N,N),BR(N,M),BI(N,M),JS(N)
        INTEGER I,J,K,IS
        REAL*8 D,P,Q,S,W

        L=1
        DO 100 K=1,N
        D=0.0
        DO 10 I=K,N
        DO 10 J=K,N
        P=AR(I,J)*AR(I,J)+AI(I,J)*AI(I,J)
        IF (P.GT.D) THEN
        D=P
        JS(K)=J
        IS=I
        END IF
10      CONTINUE
      W=D
        IF (W+1.0.EQ.1.0) THEN
        WRITE(*,20)
        L=0
        RETURN
        END IF
20    FORMAT(1X,'ERR**FAIL')
      DO 30 J=K,N
        P=AR(K,J)
        AR(K,J)=AR(IS,J)
        AR(IS,J)=P
        P=AI(K,J)
        AI(K,J)=AI(IS,J)
        AI(IS,J)=P
30      CONTINUE
      DO 35 J=1,M
        P=BR(K,J)
        BR(K,J)=BR(IS,J)
        BR(IS,J)=P
        P=BI(K,J)
        BI(K,J)=BI(IS,J)
        BI(IS,J)=P
35      CONTINUE
      DO 50 I=1,N
        P=AR(I,K)
        AR(I,K)=AR(I,JS(K))
        AR(I,JS(K))=P
        P=AI(I,K)
        AI(I,K)=AI(I,JS(K))
        AI(I,JS(K))=P
50      CONTINUE
      DO 60 J=K+1,N
        P=AR(K,J)*AR(K,K)
        Q=-AI(K,J)*AI(K,K)
        S=(AR(K,K)-AI(K,K))*(AR(K,J)+AI(K,J))
        AR(K,J)=(P-Q)/D
        AI(K,J)=(S-P-Q)/D
60      CONTINUE
      DO 65 J=1,M
        P=BR(K,J)*AR(K,K)
        Q=-BI(K,J)*AI(K,K)
        S=(AR(K,K)-AI(K,K))*(BR(K,J)+BI(K,J))
        BR(K,J)=(P-Q)/D
        BI(K,J)=(S-P-Q)/D
65      CONTINUE
      DO 90 I=1,N
        IF (I.NE.K) THEN
        DO 80 J=K+1,N
        P=AR(I,K)*AR(K,J)
        Q=AI(I,K)*AI(K,J)
        S=(AR(I,K)+AI(I,K))*(AR(K,J)+AI(K,J))
        AR(I,J)=AR(I,J)-P+Q
        AI(I,J)=AI(I,J)-S+P+Q
80      CONTINUE
      DO 85 J=1,M
        P=AR(I,K)*BR(K,J)
        Q=AI(I,K)*BI(K,J)
        S=(AR(I,K)+AI(I,K))*(BR(K,J)+BI(K,J))
        BR(I,J)=BR(I,J)-P+Q
        BI(I,J)=BI(I,J)-S+Q+P
85      CONTINUE 
        END IF
90      CONTINUE
100   CONTINUE
      DO 110 K=N,1,-1
        DO 110 J=1,M
        P=BR(K,J)
        BR(K,J)=BR(JS(K),J)
        BR(JS(K),J)=P
        P=BI(K,J)
        BI(K,J)=BI(JS(K),J)
        BI(JS(K),J)=P
110     CONTINUE
      RETURN
        END







