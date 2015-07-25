      计算六角结构的光子能带，是从J到K点的能带结构，

      

        IMPLICIT NONE
        REAL*8 A(289,289),Q(289,289),B(289),C(289)
        REAL*8 A1(289,289),Q1(289,289),B1(289),C1(289)
        REAL*8 FILING,T,TC,EPS,LENTH,RATIO
        INTEGER N,M,L,M1,N1,M2,N2,INT1,INT2,I
        COMPLEX*16 G1,G2,G3,K
        REAL*8 KX,KY,PI,TEMP,XX,J1,X,Y,EPSI1,EPSI2,slope

        OPEN(UNIT=2,FILE="sguide112.dat")
        OPEN(UNIT=3,FILE="pguide112.dat")

        EPS=0.000001
        PI=3.1415926
        
        RATIO=0.3
        FILING=2.0*PI*(RATIO**2)/SQRT(3.0)

        N=289
        XX=0.0
        EPSI1=1.0
        EPSI2=2.231**2
        slope=-1.0/sqrt(3.0)
        
        DO 1 KY=0.0, 1.0/( 2.0*sqrt(3.0) ), 0.003
        KX=0.5
        K=CMPLX(KX, KY) 
        INT1=0

        DO 2 M1=-8,8
        DO 2 N1=-8,8
        G1=CMPLX( (M1-N1)/2.0, (M1+N1)*sqrt(3.0)/2.0  )
        INT1=INT1+1
        INT2=0

        DO 3 M2=-8,8
        DO 3 N2=-8,8
        INT2=INT2+1
        G2=CMPLX( (M2-N2)/2.0, (M2+N2)*sqrt(3.0)/2.0 )
        G3=G1-G2

        X=2.0*PI*ABS(G3)*RATIO
        IF (X. EQ . XX) THEN
        Y=FILING/EPSI1+(1-FILING)/EPSI2
        ELSE
        Y=(1/EPSI1-1/EPSI2)*2.0*FILING*J1(X)/X
        END If
        A(INT1,INT2)=Y*(dreal(k+G1)*dreal(k+G2)
     $               +dimag(k+G1)*dimag(k+G2)  ) 
     
        A1(INT1,INT2)=Y*( ABS(K+G1)*ABS(K+G2) ) 
        CONTINUE
       CONTINUE

        CALL CSTRQ(A,N,Q,B,C)
        CALL CSSTQ(N,B,C,Q,EPS,L)
        
        CALL CSTRQ(A1,N,Q1,B1,C1)
        CALL CSSTQ(N,B1,C1,Q1,EPS,L)

        DO 88 M1=1,288 
        DO 88 N1=M1+1,289
        IF ( B1(M1).GT.B1(N1) ) THEN
        TEMP=B1(M1)
        B1(M1)=B1(N1)
        B1(N1)=TEMP
        END IF

        IF ( B(M1).GT.B(N1) ) THEN
        TEMP=B(M1)
        B(M1)=B(N1)
        B(N1)=TEMP
        END IF
       CONTINUE

        INT2=0
        DO 44 INT1=1,289
        IF (INT2.EQ.7) THEN
        GOTO 41
        END IF
        IF (B(INT1).GE.XX) THEN
        INT2=INT2+1
        
        WRITE(2,*) kx+ky, 2.0*SQRT(ABS(B(INT1)))/sqrt(3.0)
        END IF
       CONTINUE

       INT2=0
        DO 43 INT1=1,289
        IF (INT2.EQ.7) THEN
        GOTO 1
        END IF
        IF (B1(INT1).GE.XX) THEN
        INT2=INT2+1
        WRITE(3,*) KX+KY, 2.0*SQRT( ABS(B1(INT1)) )/sqrt(3.0)
        END IF 
       continue

        CONTINUE

        END

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

        SUBROUTINE CSTRQ(A,N,Q,B,C)
        INTEGER I,J,N
        REAL*8 A(N,N),Q(N,N),B(N),C(N)
        REAL*8 H,F,G,H2
        DO 51 I=1,N
        DO 51 J=1,N
        Q(I,J)=A(I,J)
       CONTINUE

        DO 52 I=N,2,-1
         H=0.0
        IF (I. GT .2) THEN
          DO 53 K=1,I-1
        H=H+Q(I,K)*Q(I,K)
       CONTINUE

        END IF
        IF (H+1.0. EQ .1.0) THEN
         C(I)=0.0
        IF (I. EQ .2) THEN
        C(I)=Q(I,I-1)
        END IF
        B(I)=0.0
        ELSE
         C(I)=SQRT(H)
        IF (Q(I,I-1). GT .0.0) THEN
        C(I)=-C(I)
        END IF
        H=H-Q(I,I-1)*C(I)
        Q(I,I-1)=Q(I,I-1)-C(I)
        F=0.0
        DO 50 J=1,I-1
        Q(J,I)=Q(I,J)/H
        G=0.0
        DO 54 K=1,J
        G=G+Q(J,K)*Q(I,K)
       CONTINUE
        IF (J+1. LE .I-1) THEN
        DO 55 K=J+1,I-1
        G=G+Q(K,J)*Q(I,K)
       CONTINUE
        END IF
        C(J)=G/H
        F=F+G*Q(J,I)
50      CONTINUE
        H2=F/(H+H)
        DO 56 J=1,I-1
        F=Q(I,J)
        G=C(J)-H2*F
        C(J)=G
        DO 57 K=1,J
        Q(J,K)=Q(J,K)-F*C(K)-G*Q(I,K)
57      CONTINUE
56      CONTINUE
      B(I)=H
        END IF
52      CONTINUE
        DO 58 I=1,N-1
        C(I)=C(I+1)
58      CONTINUE
        C(N)=0.0
        B(1)=0.0
        DO 59 I=1,N
         IF ((B(I). NE. 0.0).AND.(I-1.GE.1)) THEN
        DO 62 J=1,I-1
        G=0.0
        DO 60 K=1,I-1
        G=G+Q(I,K)*Q(K,J)
60      CONTINUE
        DO 61 K=1,I-1
        Q(K,J)=Q(K,J)-G*Q(K,I)
61      CONTINUE
62      CONTINUE
        END IF
        B(I)=Q(I,I)
        Q(I,I)=1.0
        IF (I-1 .GE. 1) THEN
        DO 63 J=1,I-1
         Q(I,J)=0.0
         Q(J,I)=0.0
       CONTINUE
        END IF
       CONTINUE
      RETURN
        END


        SUBROUTINE CSSTQ(N,B,C,Q,EPS,L)
        INTEGER N,L,IT
        REAL*8 B(N),C(N),Q(N,N),EPS
        REAL*8 D,H,P,R,F,E,S,G
        C(N)=0.0
        D=0.0
        F=0.0
        DO 70 J=1,N
        IT=0
        H=EPS*( ABS(B(J))+ABS(C(J)) )
        IF (H.GT.D) THEN
        D=H
        END IF
        M=J-1
       M=M+1
        IF (M.LE.N) THEN
         IF (ABS(C(M)).GT.D) GOTO 71
        END IF
         IF (M.NE.J) THEN
      IF (IT.EQ.60) THEN
      L=0
        WRITE(*,88)
       FORMAT (1X,'FAIL')
        RETURN
         END IF
        IT=IT+1
        G= B(J)
        P=(B(J+1)-G)/( 2.0*C(J) )
        R=SQRT(P*P+1.0)
        IF (P.GE.0.0) THEN
        B(J)=C(J)/(P+R)
        ELSE
        B(J)=C(J)/(P-R)
        END IF
        H=G-B(J)
        DO 72 I=J+1,N
        B(I)=B(I)-H
       CONTINUE

      F=F+H
        P=B(M)
        E=1.0
        S=0.0
        DO 74 I=M-1,J,-1
        G=E*C(I)
        H=E*P
        IF (ABS(P).GE.ABS(C(I))) THEN
        E=C(I)/P
        R=SQRT(E*E+1.0)
        C(I+1)=S*P*R
        S=E/R
        E=1.0/R
        ELSE
        E=P/C(I)
        R=SQRT(E*E+1.0)
        C(I+1)=S*C(I)*R
        S=1.0/R
        E=E/R
        END IF
        P=E*B(I)-S*G
        B(I+1)=H+S*(E*G+S*B(I))
        DO 73 K=1,N
        H=Q(K,I+1)
        Q(K,I+1)=S*Q(K,I)+E*H
        Q(K,I)=E*Q(K,I)-S*H
       CONTINUE
       CONTINUE
        C(J)=S*P
        B(J)=E*P
        IF (ABS(C(J)).GT.D) GOTO 75
        END IF
        B(J)=B(J)+F
       CONTINUE
      DO 76 I=1,N
        K=I
        P=B(I)
        IF (I+1. le .N) THEN
        J=I
       J=J+1
         IF (J.LE.N) THEN
         IF (B(J).LE.P) THEN
        K=J
        P=B(J)
        GOTO 80
          END IF
         END IF
        END IF
        IF (K.NE.I) THEN
         B(K)=B(I)
         B(I)=P
        DO 77 J=1,N
         P=Q(J,I)
         Q(J,I)=Q(J,K)
         Q(J,K)=P
       CONTINUE
      END IF
       CONTINUE
      L=1 
        RETURN
        END
