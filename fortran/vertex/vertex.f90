SUBROUTINE u4ind(u_out,rcl,l,N,TM)

  IMPLICIT NONE

  INTEGER,INTENT(in) :: l,N
  COMPLEX*16, DIMENSION(N,N), INTENT(in) :: TM   !Transformation Matrix
  DOUBLE PRECISION, DIMENSION(2*l+1,2*l+1,2*l+1,2*l+1) :: uc
  !double precision, dimension(N,N,N,N), intent(out) :: u_out
  COMPLEX*16, DIMENSION(N,N,N,N), INTENT(out) :: u_out
  DOUBLE PRECISION, DIMENSION(N,N,N,N) :: u_tmp
  DOUBLE PRECISION, DIMENSION(l+1), INTENT(in) :: rcl
  
  INTEGER :: mmax,k,k2p1,ms1,ms2,ms3,ms4,ms5,ms6,ms7,ms8
  INTEGER :: sp,ms1sig,ms2sig,ms3sig,ms4sig
  INTEGER :: xk,xm1,xm2,xm3,xm,xm4
  DOUBLE PRECISION :: cgk0,cgk1,cgk2,yor(7,7),yoi(7,7)
  COMPLEX*16 :: am1,am2,am3,am4

  !external cgk, ctormt

  !WRITE(*,*)l,N    
  mmax=2*l+1
   
  IF ((N==mmax).OR.(N==2*mmax)) THEN
     ! dimensions are fine:
     uc=0.d0

     DO k = 0, 2*l, 2
        k2p1 = k/2 + 1
        cgk0 =  cgk(l,0,k,0,l,0)
        DO ms1 = 1,mmax
           xm1 = (ms1-l-1)
           DO ms2 = 1,mmax
              xm2 = (ms2-l-1)
              DO ms3 = 1,mmax
                 xm3 = (ms3-l-1)
                 xm  = xm1 - xm3
                 DO ms4 = 1,mmax
                    IF ((ms1+ms2-ms3-ms4).NE.0) CYCLE
                    xm4 = (ms4-l-1)
                    cgk1 =  cgk(l,xm3,k,xm,l,xm1)
                    cgk2 =  cgk(l,xm2,k,xm,l,xm4)
                    uc(ms1,ms2,ms3,ms4) = uc(ms1,ms2,ms3,ms4) + rcl(k2p1)*cgk0*cgk0*cgk1*cgk2
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     u_tmp = 0.d0
     sp = (N/mmax)-1

     ! Now construct the big u matrix:
     ! expand in spins:
     DO ms1=1,mmax
        DO ms1sig =0,sp
           DO ms2=1,mmax
              DO ms2sig=0,sp
                 DO ms3 = 1,mmax
                    DO ms3sig = 0,sp
                       DO ms4 = 1,mmax
                          DO ms4sig = 0,sp
                             IF ((ms1sig==ms3sig).AND.(ms2sig==ms4sig)) THEN
                                u_tmp(ms1sig*mmax+ms1,ms2sig*mmax+ms2,ms3sig*mmax+ms3,ms4sig*mmax+ms4) = uc(ms1,ms2,ms3,ms4)
                             ENDIF
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     
     !call ctormt()
     ! Transformation:
     !write(*,*)'TEST'
     u_out = 0.d0
     DO ms1=1,N
        DO ms2=1,N
           DO ms3=1,N
              DO ms4=1,N
                 DO ms5=1,N
                    am1 = CONJG(TM(ms1,ms5))   !cmplx(yor(ms1,ms5),-yoi(ms1,ms5))
                    DO ms6=1,N
                       am2 = CONJG(TM(ms2,ms6))  !cmplx(yor(ms2,ms6),-yoi(ms2,ms6))
                       DO ms7=1,N
                          am3 = TM(ms3,ms7)      !cmplx(yor(ms3,ms7),yoi(ms3,ms7))
                          DO ms8=1,N
                             am4 = TM(ms4,ms8)    !cmplx(yor(ms4,ms8),yoi(ms4,ms8))
                             u_out(ms1,ms2,ms3,ms4) = u_out(ms1,ms2,ms3,ms4) + am1*am2*am3*am4 * u_tmp(ms5,ms6,ms7,ms8)
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
                 !if (abs(u_out(ms1,ms2,ms3,ms4))>0.0001d0) then
                 !   write(*,*)ms1,ms2,ms3,ms4, u_out(ms1,ms2,ms3,ms4)
                 !endif
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     
  ELSE

     WRITE(*,*)"N and l does not fit together: N=2l+1 or 2*(2l+1)!"

  ENDIF
  

     

CONTAINS

         
  SUBROUTINE ctormt
    IMPLICIT NONE
    DOUBLE PRECISION :: sqtwo, sq54, sq34
    
    yor=0.d0
    yoi=0.d0
    sqtwo=1.d0/SQRT(2.d0)
    sq54=SQRT(5.d0)/4d0
    sq34=SQRT(3.d0)/4d0
    IF (l.EQ.0) THEN
       yor(1,1)=1.d0
    ELSEIF (l.EQ.1) THEN
       yor(1,1)= sqtwo
       yor(1,3)=-sqtwo
       yor(2,2)=1.d0
       yoi(3,1)= sqtwo
       yoi(3,3)= sqtwo
    ELSEIF (l.EQ.2) THEN
       !yoi(1,1)= sqtwo
       !yoi(1,5)=-sqtwo
       !yoi(2,2)= sqtwo
       !yoi(2,4)= sqtwo
       !yor(3,3)=1.d0
       !yor(4,2)= sqtwo
       !yor(4,4)=-sqtwo
       !yor(5,1)= sqtwo
       !yor(5,5)= sqtwo
       ! Wien2K matrix:
       yor(3,1) = -sqtwo
       yor(3,5) =  sqtwo
       yor(5,2) =  sqtwo
       yor(5,4) =  sqtwo
       yor(1,3) =  1.d0
       yor(4,2) =  sqtwo
       yor(4,4) = -sqtwo
       yor(2,1) =  sqtwo
       yor(2,5) =  sqtwo
    ELSEIF (l.EQ.3) THEN
       yoi(1,2)=sqtwo
       yoi(1,6)=-sqtwo
       yor(2,1)=-sq54
       yor(2,3)=sq34
       yor(2,5)=-sq34
       yor(2,7)=sq54
       yoi(3,1)=-sq54
       yoi(3,3)=-sq34
       yoi(3,5)=-sq34
       yoi(3,7)=-sq54
       yor(4,4)=1.d0
       yor(5,1)=-sq34
       yor(5,3)=-sq54
       yor(5,5)=sq54
       yor(5,7)=sq34
       yoi(6,1)=sq34
       yoi(6,3)=-sq54
       yoi(6,5)=-sq54
       yoi(6,7)=sq34
       yor(7,2)=sqtwo
       yor(7,6)=sqtwo
    ENDIF
  END SUBROUTINE ctormt


  DOUBLE PRECISION FUNCTION cgk(a,al,b,be,c,ga)
    IMPLICIT NONE
    
    INTEGER :: a,al,b,be,c,ga
    
    INTEGER :: z,zmin,zmax,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13
    DOUBLE PRECISION :: fa(0:20), fac
    
    fa = (/1.d0, 1.d0, 2.d0, 6.d0, 24.d0, 12.d1, 72.d1, 504.d1,&
         4032.d1, 36288.d1, 36288.d2, 399168.d2, 4790016.d2,&
         62270208.d2, 871782912.d2, 1307674368.d3, 20922789888.d3,&
         355687428096.d3, 6402373705728.d3, 121645100408832.d3,243290200817664.d4/)
    
    
    i1=0
    i2=(a+b-c)
    i3=(a-al)
    i4=(b+be)
    i5=(c-b+al)
    i6=(c-a-be)
    zmin=MAX(i1,-i5,-i6)
    zmax=MIN(i2, i3, i4)
    cgk=0.d0
    IF (ABS(al).GT.a)   RETURN
    IF (ABS(be).GT.b)   RETURN
    IF (ABS(ga).GT.c)   RETURN
    IF ( zmin.GT.zmax )  RETURN
    IF ( (al+be).NE.ga ) RETURN
    i7=(a-b+c)
    i8=(c+b-a)
    i9=(c+b+a)
    i10=(a+al)
    i11=(b-be)
    i12=(c+ga)
    i13=(c-ga)
    DO z=zmin,zmax
       IF (MOD(z,2)==0) THEN
          fac = 1.d0
       ELSE
          fac=-1.d0
       ENDIF
       cgk=cgk+fac/(fa(z)*fa(i2-z)*fa(i3-z)*fa(i4-z)*fa(i5+z)* fa(i6+z))
    ENDDO
    cgk=cgk*SQRT(fa(i2)*fa(i7)*fa(i8)*fa(i10)*fa(i3)*fa(i4)*fa(i11)*fa(i12)*fa(i13)*(2.d0*c+1.d0)/fa(i9+1))
    
  END FUNCTION cgk

END SUBROUTINE u4ind
