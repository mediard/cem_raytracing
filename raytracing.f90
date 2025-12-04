MODULE raytracingmod
 USE ALL_DATA
 USE commandresults
 USE realloc_mod
 IMPLICIT NONE
 
CONTAINS
!888888888888888888888888888888888888888888
!888888888888888888888888888888888888888888
SUBROUTINE raytracinggMod(inputfilename)
! CALL readinput() ! This is now in indoor.f90
IMPLICIT NONE
CHARACTER(len=12), INTENT(IN) :: inputfilename
 INTEGER :: i,n,jj,npoint, pointcount,ii,j
 REAL, DIMENSION(4,R+3):: path1, path2
 INTEGER :: stat, ORD1,iii, stat1, ORD
 REAL, DIMENSION(3) :: loc_obs1
 REAL :: theta1,theta_inc, deltax, deltay
 REAL :: abs_valx, abs_valy, abs_valz
 REAL :: phase_valx, phase_valy, phase_valz
 COMPLEX :: Rpar1, Rnorm1, Tpar1, Tnorm1,Exx,Eyy,Ezz
 COMPLEX , DIMENSION (1:90,4) :: coeffs
 COMPLEX, DIMENSION(3) :: Ei1
 REAL :: pathl, Eres2nd, Eres3rd, Eres4th, EMULTI
 COMPLEX, DIMENSION(3) :: EEtemp
 COMPLEX :: EE
COMPLEX, DIMENSION(3) :: Erefl1st, Erefl1st_h, Erefl2nd_h
INTEGER :: rflpnt, rflcoord
INTEGER :: iWc, iWc_1st_h, iWc_2h
REAL :: Emean_rec0 ! eliminate
REAL :: EdandEm1, Em2nhigh, Ed0,Em1,Em1n2, Em3nhigh
CALL readinput(inputfilename)
CALL Get_data()
!CALL ReadDirectivity()
OPEN(UNIT = 10, FILE = 'Exxxxxxx.TBL', STATUS = 'UNKNOWN')
OPEN(UNIT = 20, FILE = 'Eyyyyyyy.TBL', STATUS = 'UNKNOWN')
OPEN(UNIT = 30, FILE = 'Ezzzzzzz.TBL', STATUS = 'UNKNOWN')
OPEN(UNIT = 40, FILE = 'Ettttttt.TBL', STATUS = 'UNKNOWN')
OPEN(UNIT = 50, FILE = 'Emeannnn.TBL', STATUS = 'UNKNOWN')
OPEN(UNIT = 87, FILE = 'EdandEm1.TBL', STATUS = 'UNKNOWN') ! 
OPEN(UNIT = 88, FILE = 'Em2nhigh.TBL', STATUS = 'UNKNOWN') ! 
OPEN(UNIT = 51, FILE = 'EM_MULTI.TBL', STATUS = 'UNKNOWN') 
OPEN(UNIT = 52, FILE = 'Em_2_res.TBL', STATUS = 'UNKNOWN') ! 
OPEN(UNIT = 53, FILE = 'Em_3_res.TBL', STATUS = 'UNKNOWN') ! 
OPEN(UNIT = 54, FILE = 'Em_4_res.TBL', STATUS = 'UNKNOWN') ! 
OPEN(UNIT = 70, FILE = 'Em_1_res.TBL', STATUS = 'UNKNOWN') ! 
OPEN(UNIT = 71, FILE = 'Ed_0_dir.TBL', STATUS = 'UNKNOWN') ! 
OPEN(UNIT = 72, FILE = 'Em_12res.TBL', STATUS = 'UNKNOWN') ! 
OPEN(UNIT = 73, FILE = 'Em3nhigh.TBL', STATUS = 'UNKNOWN') ! 
 WRITE(*,*) "TL = ", TL
WRITE(*,*) "R = ", R
 WRITE(*,*)All_Walls(1,:)
 CALL Get_Tree()
 WRITE(*,*) 'W =',W
 WRITE(*,*) SHAPE(Tree)
 WRITE(*,*) "TL = ", TL
 WRITE(*,"(21F4.0)") (Tree(0,i),i=0,20)

 WRITE(*,"(21F4.0)") (Tree(1,i),i=0,20)
 WRITE(*,"(21F4.0)") (Tree(2,i),i=0,20)
 WRITE(*,"(21F4.0)") (Tree(3,i),i=0,20)
 WRITE(*,"(21F4.0)") (Tree(4,i),i=0,20)
! for the moment assume area_indicator = 3
 WRITE(*,*)"TX -",TX
overx :DO i=1,xn+1
       WRITE(*,*) 'over x i=',i,' out of ',xn+1 
       WRITE(10,*) xrec(i)
       WRITE(20,*) xrec(i)
       WRITE(30,*) xrec(i)
       WRITE(40,*) xrec(i)
       WRITE(50,*) xrec(i)
       WRITE(51,*) xrec(i)
       WRITE(52,*) xrec(i)
       WRITE(53,*) xrec(i)
       WRITE(54,*) xrec(i)
       WRITE(87,*) xrec(i)
       WRITE(88,*) xrec(i)
       WRITE(70,*) xrec(i)
       WRITE(71,*) xrec(i)
       WRITE(72,*) xrec(i)
       WRITE(73,*) xrec(i)
       overy: DO j=1,yn+1
        loc_obs1=(/ xrec(i), yrec(j), zrec(1) /)
       Ei1 = (/ (0.,0.) , (0.,0.) , (0.,0.) /) 
       Emean_rec0 = 0.
       Eres2nd = 0.
       Eres3rd = 0.
       Eres4th = 0.
       EMULTI = 0.

       iii = 1
       ii = 0
                         CALL Get_one_path(ii, loc_obs1, path2, stat, ORD) ! change to if stat==0 then OK
             pathstatus:    IF (stat == 0) THEN
                 CALL Path_length(path2, pathl)
                 Ei1 = Path_field (path2,pathl)
  !               Ert_dir_refl1st = Ei1
                 Emean_rec0 = ABS(Ei1(1))**2+ABS(Ei1(2))**2+ABS(Ei1(3))**2
                 IF (ORD==0) THEN
                         EdandEm1=Emean_rec0
                         Ed0 = EdandEm1
                         Em1n2 = 0.
                         Em3nhigh=0.
                         Em2nhigh=0.
                 ELSE
                         STOP 'error in order of direct field'
                 END IF
                 iWc = 0
                 iWc_1st_h = 0
                 iWc_2h = 0
         END IF pathstatus

              fl: DO 
                     ii = ii+1
                     ! WRITE(*,*) 'ii = ', ii
                     IF (ii>TL) EXIT fl
                     CALL Get_one_path(ii,loc_obs1,path2, stat1 ,ORD)
                     pathstatus2: IF (stat1 == 0) THEN
                     IF (path2(1,2)==0. .AND. path2(1,3)==0.) THEN
                             WRITE(*,*) 'exiting fl from 0 0 '
                             EXIT fl
                     END IF
                     iii = iii+1
                     CALL Path_length(path2, pathl)
                     EEtemp=Path_field( path2 , pathl )
                     Emean_rec0 = Emean_rec0+ &
                            ABS(EEtemp(1))**2+ABS(EEtemp(2))**2+ABS(EEtemp(3))**2
                     Ei1 = Ei1 + EEtemp
                    IF (ORD>=2) THEN
                            Em2nhigh = Em2nhigh + &
                            ABS(EEtemp(1))**2+ABS(EEtemp(2))**2+ABS(EEtemp(3))**2
                    END IF
                     IF (ORD==4) THEN
                             Eres4th = Eres4th + &
                            ABS(EEtemp(1))**2+ABS(EEtemp(2))**2+ABS(EEtemp(3))**2
                    END IF
                    IF (ORD==3) THEN
                            Eres3rd = Eres3rd + &
                            ABS(EEtemp(1))**2+ABS(EEtemp(2))**2+ABS(EEtemp(3))**2
                    END IF
                    IF (ORD>=3) THEN
                            Em3nhigh = Em3nhigh + &
                            ABS(EEtemp(1))**2+ABS(EEtemp(2))**2+ABS(EEtemp(3))**2
                    END IF
                    IF (ORD==2) THEN
                            Eres2nd = Eres2nd + &
                            ABS(EEtemp(1))**2+ABS(EEtemp(2))**2+ABS(EEtemp(3))**2
                            Em1n2 = Em1n2 + &
                            ABS(EEtemp(1))**2+ABS(EEtemp(2))**2+ABS(EEtemp(3))**2
                    END IF
                    IF (ORD>=1) THEN
                            EMULTI = EMULTI + &
                            ABS(EEtemp(1))**2+ABS(EEtemp(2))**2+ABS(EEtemp(3))**2
                    END IF
                    IF (ORD==1) THEN
                            EdandEm1 = EdandEm1 + &
                            ABS(EEtemp(1))**2+ABS(EEtemp(2))**2+ABS(EEtemp(3))**2
                            Em1 = Em1 + &
                            ABS(EEtemp(1))**2+ABS(EEtemp(2))**2+ABS(EEtemp(3))**2
                            Em1n2 = Em1n2 + &
                            ABS(EEtemp(1))**2+ABS(EEtemp(2))**2+ABS(EEtemp(3))**2
                    END IF

             END IF pathstatus2
              END DO fl
              n = ii-1
              Exx=Ei1(1)
              Eyy=Ei1(2)
              Ezz=Ei1(3)
              ! for x
              abs_valx = ABS(Exx)
              IF (REAL(Exx)>=0) THEN
                     phase_valx = atan(AIMAG(Exx)/REAL(Exx)) *180./pi
                     IF (phase_valx<0.) phase_valx=phase_valx+360.
              ELSEIF (REAL(Exx)<0) THEN
                     phase_valx = atan(AIMAG(Exx)/REAL(Exx))
                     IF (phase_valx>=0.) phase_valx=phase_valx+180.
                     IF (phase_valx<0.) phase_valx=180.-phase_valx
              ENDIF
              WRITE(10,*) yrec(j), abs_valx/SQRT(2.)
              abs_valy = ABS(Eyy)
              IF (REAL(Eyy)>=0) THEN
                     phase_valy = atan(AIMAG(Eyy)/REAL(Eyy)) *180./pi
                     IF (phase_valy<0.) phase_valy=phase_valy+360.
              ELSEIF (REAL(Eyy)<0) THEN
                     phase_valy = atan(AIMAG(Eyy)/REAL(Eyy))
                     IF (phase_valy>=0.) phase_valy=phase_valy+180.
                     IF (phase_valy<0.) phase_valy=180.-phase_valy
              ENDIF
              WRITE(20,*) yrec(j), abs_valy/SQRT(2.)
              abs_valz = ABS(Ezz)
              IF (REAL(Ezz)>=0) THEN
                     phase_valz = atan(AIMAG(Ezz)/REAL(Ezz)) *180./pi
                     IF (phase_valz<0.) phase_valz=phase_valz+360.
              ELSEIF (REAL(Ezz)<0) THEN
                     phase_valz = atan(AIMAG(Ezz)/REAL(Ezz))
                     IF (phase_valz>=0.) phase_valz=phase_valz+180.
                     IF (phase_valz<0.) phase_valz=180.-phase_valz
              ENDIF
              WRITE(30,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(abs_valx**2+abs_valy**2+abs_valz**2)
              WRITE(40,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(Emean_rec0)
              WRITE(50,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(EMULTI)
              WRITE(51,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(Eres2nd)
              WRITE(52,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(Eres3rd)
              WRITE(53,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(Eres4th)
              WRITE(54,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(EdandEm1)
              WRITE(87,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(Em2nhigh)
              WRITE(88,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(Em1)
              WRITE(70,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(Ed0)
              WRITE(71,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(Em1n2)
              WRITE(72,*) yrec(j), abs_valz/SQRT(2.)
              abs_valz = SQRT(Em3nhigh)
              WRITE(73,*) yrec(j), abs_valz/SQRT(2.)
       ENDDO overy
ENDDO overx

 CLOSE(10)
 CLOSE(20)
 CLOSE(30)
 CLOSE(40)
 CLOSE(50)
 CLOSE(51)
 CLOSE(52)
 CLOSE(53)
 CLOSE(54)
 CLOSE(87)
 CLOSE(88)
 CLOSE(70)
 CLOSE(71)
 CLOSE(72)
 CLOSE(73)
 WRITE(*,*) 'W = ',W
 WRITE(*,*) SHAPE(Walls)
 DO i=1,W
         WRITE(*,*) (Walls(i,j),j=1,7)
 END DO
 END SUBROUTINE raytracinggMod
 END MODULE raytracingmod
 !88888888888888888888888888888
