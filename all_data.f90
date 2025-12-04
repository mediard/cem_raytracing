MODULE ALL_DATA
REAL, PARAMETER :: pi = 3.14159265
REAL, PARAMETER :: mu0 = 4*pi*1E-7
REAL, PARAMETER :: epsilon0 = 8.854187817620E-12
REAL, PARAMETER :: c0 = 299792458
INTEGER :: W !number of walls        should be parameter
INTEGER :: n_wt  !number of wall types    parameter
INTEGER :: n_mat  ! number of materials   parameter
REAL, ALLOCATABLE,DIMENSION (:,:) :: Walls              !should be allocatable
REAL, ALLOCATABLE,DIMENSION (:,:) :: All_Walls         !should be allocatable

REAL, ALLOCATABLE,DIMENSION (:,:) :: WTypes         !should be allocatable
REAL, ALLOCATABLE,DIMENSION (:,:) :: Materials      !should be allocatable
REAL, DIMENSION (6) :: TX ! (x,y,z,power, freq, antenna_type)
INTEGER, PARAMETER :: R=6 ! maximum number of reflections
REAL, ALLOCATABLE, DIMENSION(:,:) :: Tree
INTEGER :: loctr, TL
REAL, DIMENSION(3) :: loc_tx
COMPLEX, PARAMETER :: cj = (0. , 1.)
!**********************
  REAL,DIMENSION(8,3) :: material_mat
  REAL,DIMENSION(10,15):: wall_type_mat
  CHARACTER(len=10),DIMENSION(8) :: material_names,wall_type_names
  CHARACTER(len=10),DIMENSION(200) :: point_names
  INTEGER, DIMENSION(100,2)::define_wall_mat !refers to number of points in "point_mat"
  INTEGER, DIMENSION(100) :: define_wall_types !number "wall_type_names"
  REAL,DIMENSION(200,3) :: point_mat
   CHARACTER(len=120) :: line1
   INTEGER :: pnt_num,wlltp_num,wll_def_num
   REAL, DIMENSION(0:180) :: Directivity_matrix
!**********************
 REAL,ALLOCATABLE,DIMENSION(:) :: xrec,yrec,zrec ! in meters
 INTEGER :: area_indicator ! if 1 (y,z), 2(x,z), 3(x,y)
 INTEGER :: xn,yn,zn
 !REAL,ALLOCATABLE,DIMENSION(:,:) :: Ex_abs_rec ! eliminate
 !REAL,ALLOCATABLE,DIMENSION(:,:) :: Ex_phase_rec ! eliminate
 !REAL,ALLOCATABLE,DIMENSION(:,:) :: Ey_abs_rec ! eliminate
 !REAL,ALLOCATABLE,DIMENSION(:,:) :: Ey_phase_rec ! eliminate
 !REAL,ALLOCATABLE,DIMENSION(:,:) :: Ez_abs_rec ! eliminate
 !REAL,ALLOCATABLE,DIMENSION(:,:) :: Ez_phase_rec ! eliminate
 !REAL,ALLOCATABLE,DIMENSION(:,:) :: Et_abs_rec,ESab_rec ! eliminate
 REAL :: delta_incr !in cm
!REAL, DIMENSION(:,:), ALLOCATABLE :: EmHSab,EdHSab !eliminate
!REAL, DIMENSION(:,:), ALLOCATABLE :: EmSab,EdSab ! eliminate
!COMPLEX, DIMENSION(:,:), ALLOCATABLE :: Ert_dir_refl1st ! eliminate
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_image( source_loc, wallinfo, n_vec, image_loc)
IMPLICIT NONE
REAL, DIMENSION(3), INTENT(IN) :: source_loc
REAL, DIMENSION(14), INTENT(IN) :: wallinfo
REAL, INTENT(OUT) :: n_vec
REAL, DIMENSION(3), INTENT(OUT) :: image_loc


IF (wallinfo(13)==1.) THEN
       image_loc(1) = wallinfo(1) - ( (source_loc(1))-(wallinfo(1)) )
       image_loc(2) = source_loc(2)
       image_loc(3) = source_loc(3)
       n_vec = 1. * ((source_loc(1))-(wallinfo(1)))/ABS((source_loc(1))-(wallinfo(1)))
ELSE IF (wallinfo(13)==2.) THEN
       image_loc(2) = wallinfo(2) - ( (source_loc(2))-(wallinfo(2)) )
       image_loc(1) = source_loc(1)
       image_loc(3) = source_loc(3)
       n_vec = 2. * ((source_loc(2))-(wallinfo(2)))/ABS((source_loc(2))-(wallinfo(2)))
ELSE IF (wallinfo(13)==3.) THEN
       image_loc(3) = wallinfo(3) - ( (source_loc(3))-(wallinfo(3)) )
       image_loc(1) = source_loc(1)
       image_loc(2) = source_loc(2)
       n_vec = 3. * ((source_loc(3))-(wallinfo(3)))/ABS((source_loc(3))-(wallinfo(3)))
END IF
RETURN
END SUBROUTINE Get_image
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_Data()
!USE ALL_DATA
! use the shared module; Walls, loc_tx, power(mW), Frequency(GHz), antenna_type(1)
INTEGER :: i
REAL :: temp



DO i=1, n_mat
   temp = Materials(i,2)
   Materials(i,2) = Materials(i,3)
   Materials(i,3) = temp
END DO

 CALL Process_Walls()

END SUBROUTINE Get_Data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Process_Walls()
!USE ALL_DATA
REAL :: xs,ys,zs,xe,ye,ze,wt,norm, cord
INTEGER i
ALLOCATE(All_Walls(W,14))! :: All_Walls         !should be allocatable
DO i=1,W
    xs = Walls(i,1)
    ys = Walls(i,2)
    zs = Walls(i,3)
    xe = Walls(i,4)
    ye = Walls(i,5)
    ze = Walls(i,6)
    IF (xs < xe) THEN
        xmin = xs
        xmax = xe
    ELSE IF (xs > xe) THEN
        xmin = xe
        xmax = xs
    ELSE
        xmin = xe
        xmax = xe
        norm = 1.
        cord = xe
    END IF
    IF (ys < ye) THEN
        ymin = ys
        ymax = ye
    ELSE IF (ys > ye) THEN
        ymin = ye
        ymax = ys
    ELSE
        ymin = ye
        ymax = ye
        norm = 2.
        cord = ymin
    END IF

    IF (zs < ze) THEN
        zmin = zs
        zmax = ze
    ELSE IF (zs > ze) THEN
        zmin = ze
        zmax = zs
    ELSE
        zmin = ze
        zmax = ze
        norm = 3.
        cord = ze
    END IF
All_Walls(i, 1:6) = (/xs,ys,ze,xe,ye,ze/)
All_Walls(i, 7:14) = (/xmin,ymin,zmin,xmax,ymax,zmax,norm,cord/)
END DO

END SUBROUTINE Process_Walls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_Tree()
!USE ALL_DATA
INTEGER :: wll_count, i, j, r_ord, &
           icurrent, scurrent,  last
REAL :: n_vec
REAL, DIMENSION (3) :: loc_img
! TL = Sigma_i_R W**i
TL=0
DO i=1,R
    DO j=1,W**i
        TL = TL+1
    END DO
END DO
! WRITE(*,*) "inside TL = ", TL
ALLOCATE(Tree(0:4, 0:TL))
Tree(1:3,0) = loc_tx
! loc_img
! s = -+1, -+2, -+3  ==>= -+x, -+y, -+z
last = 0
IF (R >=1) THEN
    DO wll_count = 1, W
        CALL Get_image( Tree(1:3,0), All_Walls(wll_count,:), n_vec, loc_img)
        Tree(0, wll_count) = wll_count
        Tree(1:3, wll_count) = loc_img
        Tree(4, wll_count) = n_vec
    ENDDO
END IF

IF (R>=2) THEN
    DO r_ord = 2,R
        icurrent = last + W**(r_ord-1)
            DO scurrent = (last+1), (last+W**(r_ord-1))
                DO wll_count = 1, W
                    CALL Get_image(Tree(1:3,scurrent),All_Walls(wll_count,:),n_vec, loc_img)
                    icurrent = icurrent + 1
                    Tree(0, icurrent) = wll_count
                    Tree(1:3, icurrent) = loc_img
                    Tree(4, icurrent) = n_vec
                END DO
            END DO
            last = last + W**(r_ord-1)
        ENDDO
    END IF

END SUBROUTINE Get_Tree
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE findparent(childloctr,parent_loctr,statg)
!USE ALL_DATA
IMPLICIT NONE

INTEGER :: i,wt,diff,prev,parent,temp
INTEGER, INTENT(IN) :: childloctr
INTEGER, INTENT(OUT) :: parent_loctr, statg
! childloctr : column number in Tree: INTEGER
! W: number of walls and it shared

IF (childloctr <= W) THEN
    !WRITE(*,*)'childloctr error'
    ! get out
    parent_loctr = 0
ELSE IF (childloctr > W) THEN
    i = 1
    wt = W**i
    DO
        i = i+1
        wt = wt + W**i
        IF (childloctr <= wt) EXIT
    END DO
    diff = childloctr - (wt-W**i)
    prev = FLOOR( 1.*REAL(diff-1)/REAL(W) ) + 1
    parent = wt -W**i - W**(i-1) + prev
    parent_loctr = parent
END IF
    IF ( Tree(0,childloctr)==Tree(0,parent_loctr) ) THEN
       statg = 1
    ELSE 
       statg = 0
END IF
    RETURN


END SUBROUTINE findparent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Find_inc_pt(loc_obs, temploctr, loc_inc, stat)
!USE ALL_DATA
IMPLICIT NONE
REAL, DIMENSION(3),INTENT(IN):: loc_obs
INTEGER, INTENT (IN):: temploctr 
REAL, DIMENSION(3), INTENT(OUT) :: loc_inc
INTEGER, INTENT(OUT) :: stat
REAL, DIMENSION(3) :: loc_img
REAL :: xi,yi,zi,xo,yo,zo,w_vec,x_inc,y_inc,z_inc,t
INTEGER:: w_num
REAL, DIMENSION(6) :: w_min_max
INTEGER :: myStat
REAL :: zzmin,zzmax,xxmin,xxmax,yymin,yymax
! Declarations
myStat = 0
loc_img = Tree(1:3, temploctr)
xi = loc_img(1)
yi = loc_img(2)
zi = loc_img(3)
w_num = INT(Tree(0,temploctr))
w_vec = Tree(4,temploctr)
xo = loc_obs(1)
yo = loc_obs(2)
zo = loc_obs(3)
! (x_inc, y_inc, z_inc) = (xo, yo, zo) + t * ( xi-xo , yi-yo, zi-zo )
IF (ABS(w_vec) == 1.) THEN
    x_inc = All_Walls(w_num, 1)
    IF (xi/=xo) THEN
    t = (x_inc - xo ) / (xi - xo)
    y_inc = yo + t * (yi - yo)
    z_inc = zo + t * (zi - zo)
    ELSE
            y_inc = (yi+yo)/2.
            z_inc = (zi+zo)/2.
    END IF
    !IF (xi == xo) THEN
    !        myStat=1
    !END IF
ELSE IF (ABS(w_vec) == 2.) THEN
    y_inc = All_Walls(w_num, 2)
    IF (yi/=yo) THEN !STOP 'Error yi=yo'
    t = (y_inc - yo ) / (yi - yo)
    x_inc = xo + t * (xi - xo)
    z_inc = zo + t * (zi - zo)
    ELSE
            x_inc = (xi+xo)/2.
            z_inc = (zi+zo)/2.
    END IF
    !IF (yi==yo) THEN
    !        myStat=1
    !END IF
ELSE IF (ABS(w_vec) == 3.) THEN
    z_inc = All_Walls(w_num, 3)
    IF (zi/=zo) THEN !STOP 'Error yi=yo'
    t = (z_inc - zo ) / (zi - zo)
    y_inc = yo + t * (yi - yo)
    x_inc = xo + t * (xi - xo)
    ELSE
            y_inc = (yi+yo)/2.
            x_inc = (xi+xo)/2.
    END IF
    !IF (zi==zo) THEN
    !        myStat=1
    !END IF
END IF
w_min_max = All_Walls(w_num, 7:12)
xxmin = MINVAL(All_Walls(:,7))
yymin = MINVAL(All_Walls(:,8))
zzmin = MINVAL(All_Walls(:,9))
xxmax = MINVAL(All_Walls(:,10))
yymax = MINVAL(All_Walls(:,11))
zzmax = MINVAL(All_Walls(:,12))
!
IF (myStat ==1) THEN
        stat = 1
END IF
!
IF (w_min_max(1)==w_min_max(4)) THEN
      IF (y_inc<=w_min_max(2) .OR. z_inc<=w_min_max(3) .OR. y_inc>=w_min_max(5) .OR. z_inc>=w_min_max(6)) THEN
             stat=1
     ELSE
        stat = 0
        loc_inc = (/ x_inc, y_inc, z_inc /)
      END IF
ELSE IF (w_min_max(2)==w_min_max(5)) THEN
      IF (x_inc<=w_min_max(1) .OR. z_inc<=w_min_max(3) .OR. x_inc>=w_min_max(4) .OR. z_inc>=w_min_max(6)) THEN
             stat=1
     ELSE
        stat = 0
        loc_inc = (/ x_inc, y_inc, z_inc /)
      END IF
ELSE IF (w_min_max(3)==w_min_max(6)) THEN
      IF (y_inc<=w_min_max(2) .OR. x_inc<=w_min_max(1) .OR. y_inc>=w_min_max(5) .OR. x_inc>=w_min_max(4)) THEN
              stat=1
     ELSE
        stat = 0
        loc_inc = (/ x_inc, y_inc, z_inc /)
      END IF
!END IF
!IF ( myStat==1 .OR. (x_inc<w_min_max(1)) .OR. (y_inc<w_min_max(2)) .OR. (z_inc<w_min_max(3)) .OR. &
!      (x_inc>w_min_max(4)) .OR. (y_inc>w_min_max(5)) .OR. (z_inc>w_min_max(6)) ) THEN 
      ! the >= sign instead of > helps eliminate the possibility of having inc.
      ! points on the cross section of two planes. However this will cause not
      ! accounting for any reflections. So we have to fix this somewhere else. 
!        stat = 1
END IF
!write(*,*)'loc_obs= ',loc_obs
!write(*,*)'loc_img= ',loc_img
!write(*,*)'wall = ', w_num
!write(*,*)
!write(*,*)'incdient = ',loc_inc
!write(*,*)
RETURN

END SUBROUTINE Find_inc_pt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Find_inc_pt0(loc_obs, temploctr, loc_inc, stat)
!USE ALL_DATA
IMPLICIT NONE
REAL, DIMENSION(3),INTENT(IN):: loc_obs
INTEGER, INTENT (IN):: temploctr 
REAL, DIMENSION(3), INTENT(OUT) :: loc_inc
INTEGER, INTENT(OUT) :: stat
REAL, DIMENSION(3) :: loc_img
REAL :: xi,yi,zi,xo,yo,zo,w_vec,x_inc,y_inc,z_inc,t
REAL :: zzmin,zzmax,xxmin,xxmax,yymin,yymax
INTEGER:: w_num
REAL, DIMENSION(6) :: w_min_max
INTEGER :: myStat
! Declarations
myStat = 0
loc_img = Tree(1:3, temploctr)
xi = loc_img(1)
yi = loc_img(2)
zi = loc_img(3)
w_num = INT(Tree(0,temploctr))
w_vec = Tree(4,temploctr)
xo = loc_obs(1)
yo = loc_obs(2)
zo = loc_obs(3)
! (x_inc, y_inc, z_inc) = (xo, yo, zo) + t * ( xi-xo , yi-yo, zi-zo )
IF (ABS(w_vec) == 1.) THEN
    x_inc = All_Walls(w_num, 1)
    t = (x_inc - xo ) / (xi - xo)
    y_inc = yo + t * (yi - yo)
    z_inc = zo + t * (zi - zo)
    !IF (xi == xo) THEN
    !        myStat=1
    !END IF
ELSE IF (ABS(w_vec) == 2.) THEN
    y_inc = All_Walls(w_num, 2)
    t = (y_inc - yo ) / (yi - yo)
    x_inc = xo + t * (xi - xo)
    z_inc = zo + t * (zi - zo)
    !IF (yi==yo) THEN
    !        myStat=1
    !END IF
ELSE IF (ABS(w_vec) == 3.) THEN
    z_inc = All_Walls(w_num, 3)
    t = (z_inc - zo ) / (zi - zo)
    y_inc = yo + t * (yi - yo)
    x_inc = xo + t * (xi - xo)
    !IF (zi==zo) THEN
    !        myStat=1
    !END IF
END IF
w_min_max = All_Walls(w_num, 7:12)
xxmin = MINVAL(All_Walls(:,7))
yymin = MINVAL(All_Walls(:,8))
zzmin = MINVAL(All_Walls(:,9))
xxmax = MINVAL(All_Walls(:,10))
yymax = MINVAL(All_Walls(:,11))
zzmax = MINVAL(All_Walls(:,12))

!
IF (myStat ==1) THEN
        stat = 1
END IF
!
IF (w_min_max(1)==w_min_max(4)) THEN
      IF (y_inc<=w_min_max(2) .OR. z_inc<=w_min_max(3) .OR. y_inc>=w_min_max(5) .OR. z_inc>=w_min_max(6)) THEN
             stat=1
     ELSE
        stat = 0
        loc_inc = (/ x_inc, y_inc, z_inc /)
      END IF
ELSE IF (w_min_max(2)==w_min_max(5)) THEN
      IF (x_inc<=w_min_max(1) .OR. z_inc<=w_min_max(3) .OR. x_inc>=w_min_max(4) .OR. z_inc>=w_min_max(6)) THEN
             stat=1
     ELSE
        stat = 0
        loc_inc = (/ x_inc, y_inc, z_inc /)
      END IF
ELSE IF (w_min_max(3)==w_min_max(6)) THEN
      IF (y_inc<=w_min_max(2) .OR. x_inc<=w_min_max(1) .OR. y_inc>=w_min_max(5) .OR. x_inc>=w_min_max(4)) THEN
              stat=1
     ELSE
        stat = 0
        loc_inc = (/ x_inc, y_inc, z_inc /)
      END IF
!END IF
!IF ( myStat==1 .OR. (x_inc<w_min_max(1)) .OR. (y_inc<w_min_max(2)) .OR. (z_inc<w_min_max(3)) .OR. &
!      (x_inc>w_min_max(4)) .OR. (y_inc>w_min_max(5)) .OR. (z_inc>w_min_max(6)) ) THEN 
      ! the >= sign instead of > helps eliminate the possibility of having inc.
      ! points on the cross section of two planes. However this will cause not
      ! accounting for any reflections. So we have to fix this somewhere else. 
!        stat = 1
END IF
IF (ABS(x_inc-xxmin)<1E-4 .OR. ABS(x_inc-xxmax)<1E-4 .OR. ABS(y_inc-yymin)<1E-4 .OR. &
    ABS(y_inc-yymax)<1E-4 .OR.  ABS(z_inc-zzmin)<1E-4 .OR. ABS(z_inc-zzmax)<1E-4) THEN
    stat = 1
END IF

!write(*,*)'loc_obs= ',loc_obs
!write(*,*)'loc_img= ',loc_img
!write(*,*)'wall = ', w_num
!write(*,*)
!write(*,*)'incdient = ',loc_inc
!write(*,*)
RETURN

END SUBROUTINE Find_inc_pt0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_one_path(temploctr, loc_obs, path, stat, ORD)
!USE ALL_DATA
IMPLICIT NONE
INTEGER, INTENT(IN) :: temploctr
INTEGER, INTENT(OUT) :: ORD,stat
REAL, DIMENSION(3), INTENT(IN) :: loc_obs
REAL, DIMENSION(4,R+3), INTENT(OUT) :: path
REAL, DIMENSION(4,R+3) :: path4
REAL, DIMENSION(4,R+1) :: temp_path
! Declarations
! CAREFUL with stat ... change .. now 0 means no error
INTEGER :: n, statinc,iii
REAL, DIMENSION(3) :: current_obs, current_img, loc_inc
INTEGER :: current_tr, tempcur, statg1,jj,ii
stat=0
path4=0.
temp_path=0.
n = 1
current_obs = loc_obs
IF (temploctr /= 0) THEN ! give erro msg for others
    current_img = Tree(1:3, temploctr)
    current_tr = temploctr
    wl: DO
        CALL Find_inc_pt(current_obs, current_tr, loc_inc, statinc)
        ! The condition below is for the following
        ! 1) Rejects image of image with respect to the same wall.
        ! 2) Rejects image of image with respect to wall adjacent to it. pnt of
        ! inc being on the intersection. wall are parallel.
        ! 3) Reflection from edges and corners.
        !IF ((current_obs(1)==loc_inc(1)) .AND. (current_obs(2)==loc_inc(2)) .AND. (current_obs(3)==loc_inc(3)) ) THEN
        IF (ABS(current_obs(1)-loc_inc(1)) + ABS(current_obs(2)-loc_inc(2))&
        + ABS(current_obs(3)-loc_inc(3)) < 1E-6) THEN
                !WRITE(*,*) 'Find Identical Point ERROR'
                !READ(*,*)
                stat = 1 
                EXIT wl
        ELSEIF (statinc == 0) THEN ! careful with 0
            temp_path(1,n) = current_tr
            temp_path(2:4,n) = loc_inc
            n = n+1
            current_obs = loc_inc
            tempcur = current_tr
            CALL findparent(tempcur,current_tr,statg1)
            IF (statg1==1) THEN
                  stat = 1
                  !WRITE(*,*) 'FIND PARENT ERROR'
                  !READ(*,*)
                  EXIT wl
            END IF
        ELSE IF (statinc == 1) THEN
            stat = 1
            !WRITE(*,*) 'Find INC. Point ERROR'
            !READ(*,*) 
            EXIT wl
        END IF
        IF (current_tr ==0) THEN  
            temp_path(1,n) = 0.!!!!!!!!!!!!!!!!!!!!!!
            temp_path(2:4,n) = loc_tx!!!!!!!!!!!!!!!!!!
            EXIT wl
        END IF
    END DO wl
    IF (stat==0) THEN
        DO ii=1,4
            DO jj=1,(R+1)
                path4(ii,2+jj)=temp_path(ii,jj)
            END DO
        END DO
        !path(1:4,3:(R+3)) = temp_path
        stat = 0
    END IF
    ORD = n-1
ELSE IF (temploctr==0) THEN
    !path4(1:4,3) = (/0., loc_tx(1),loc_tx(2),loc_tx(3) /) ! check if it's possilbe
    path4(1,3) =0.
    path4(2,3) = loc_tx(1)
    path4(3,3) = loc_tx(2)
    path4(4,3) = loc_tx(3)
    ORD = n-1
!    WRITE(*,*) 'zero ORD = ', ORD
!    READ(*,*)
END IF

IF (stat==0) THEN
path4(1:4,1) = REAL(ORD)

path4(1,2) = 0.
path4(2:4,2) = loc_obs

path=path4

END IF
                 !   IF (stat==0) THEN
                 !          WRITE(*,*) "path is "
                 !          WRITE(*,*) (path(1,jj),jj=1,(R+3))
                 !          WRITE(*,*) (path(2,jj),jj=1,(R+3))
                 !          WRITE(*,*) (path(3,jj),jj=1,(R+3))
                 !          WRITE(*,*) (path(4,jj),jj=1,(R+3))
                 !          WRITE(*,*) "statinc = ", statinc
                 !          WRITE(*,*) "statg1 = ", statg1
                 !          WRITE(*,*) "ORD = ", path(1,1)
                 !          READ(*,*)
                 !  END IF
!DO iii=1,(R+3)
!IF (ORD == 0) THEN
!END DO
RETURN
END SUBROUTINE Get_one_path
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE Get_all_paths(loc_obs, All_Paths)
! !USE ALL_DATA
!IMPLICIT NONE
!REAL, DIMENSION(3), INTENT(IN) :: loc_obs
!REAL, DIMENSION(TL+1,R+3,4),INTENT(OUT) :: All_Paths
!REAL, DIMENSION(TL+1,R+3,4) :: temp
!INTEGER :: i, stat, ORD,n,jj,jjj
!REAL, DIMENSION(4,R+3) :: path2
!
!! All_Paths : (length_tr) X (R+3) X (1:4)
!n = 0
!All_Paths = 0.
!!OPEN(UNIT = 14, FILE = 'ALLPATHS.dat', STATUS = 'UNKNOWN')
!
!DO i=0, TL
!    
!    CALL Get_one_path(i, loc_obs, path2, stat, ORD) ! change to if stat==0 then OK
!  
!    IF (stat == 0) THEN
!        n=n+1
!        temp(n,:,1:4) = TRANSPOSE(path2) ! CAREFUL
!        
!    END IF
!
!END DO
!
!All_Paths(1:n,:,:) = temp
!
!!WRITE(*,*) "ALL_Paths from inside:",All_Paths
!
!RETURN
!END SUBROUTINE Get_all_paths
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_angle(loc_obs, r_dist, theta, phi)
!USE ALL_DATA
IMPLICIT NONE
REAL, DIMENSION(3), INTENT(IN) :: loc_obs
!REAL, DIMENSION(3), INTENT(IN) :: loc_tx
REAL, INTENT(OUT) :: r_dist, theta, phi
REAL :: x,y,z
x = loc_obs(1) - loc_tx(1)
y = loc_obs(2) - loc_tx(2)
z = loc_obs(3) - loc_tx(3)
r_dist = SQRT(x**2 + y**2 + z**2)
IF (x/=0. .OR. y/=0.) THEN
       IF (y>=0.) THEN
              phi = ACOS(x/SQRT(x**2+y**2))
       ELSE
              phi = 2.*pi - ACOS(x/SQRT(x**2+y**2))
       END IF
ELSE
        phi = 0. ! put anything.
END IF
IF (z==0. .AND. r_dist==0.) THEN
        theta = 1. ! this is error. on the antenna itself
ELSE
        theta = ACOS(z/r_dist)
END IF
RETURN
END SUBROUTINE Get_angle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION Directivity(theta, phi, ant_type)
!USE ALL_DATA
IMPLICIT NONE
REAL, INTENT(IN) :: theta, phi
INTEGER, INTENT(IN) :: ant_type
REAL :: Directivity
IF (ant_type /= 1) THEN
    WRITE(*,*) 'antenna type error'
ELSE IF (theta/=0.) THEN
    Directivity = COS(pi/2.*COS(theta))/SIN(theta)
ELSE
        Directivity = 0.
END IF
END FUNCTION Directivity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadDirectivity()
IMPLICIT NONE
INTEGER :: i
REAL :: d_theta,d_phi, d_d
CHARACTER(len=200) :: d_line 
OPEN(UNIT=201, FILE='directivity.txt', STATUS='UNKNOWN')
        READ(201,'(A)') d_line
        READ(201,'(A)') d_line
DO i=0,180
        READ(201,'(A)') d_line
        READ(d_line,*) d_theta, d_phi, d_d
        Directivity_matrix(i) = d_d
END DO
Directivity_matrix = Directivity_matrix/MAXVAL(Directivity_matrix)
DO i=0,180
        WRITE(*,*) Directivity_matrix(i)
END DO
READ(*,*)
END SUBROUTINE ReadDirectivity
!88888888888888888888888888888888888888888888888888888
!88888888888888888888888888888888888888888888888888888
FUNCTION DirectivityData(theta, phi, ant_type)
!USE ALL_DATA
IMPLICIT NONE
REAL, INTENT(IN) :: theta, phi
INTEGER, INTENT(IN) :: ant_type
REAL :: DirectivityData
DirectivityData = Directivity_matrix(FLOOR(theta*180./pi))
END FUNCTION DirectivityData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE R_T_coeffs (theta, freq, w_num, n_lay, Rpar, Rnorm, Tpar, Tnorm)
!USE ALL_DATA
IMPLICIT NONE
REAL, INTENT(IN) :: theta, freq
INTEGER, INTENT(IN) :: w_num,  n_lay
COMPLEX, INTENT(OUT) :: Rpar, Rnorm, Tpar, Tnorm

REAL :: omega

REAL, DIMENSION(0:n_lay+1) :: loss_tans, murs, epsilonrs, sigmas
COMPLEX, DIMENSION(0:n_lay+1) :: Gammas, thetas, A, B, C, D, Y, Z
COMPLEX, DIMENSION(0:n_lay) :: psi
COMPLEX :: Gamma0,temp, ptemp, mtemp
INTEGER :: i, wt, j, statc

wt = INT(Walls(w_num, 7))
!WRITE(*,*) 'wt = ' ,wt

omega = 2*pi*freq
loss_tans(0) = 0.
! REAL --> COMPLEX
murs(0) = 1.
epsilonrs(0) = 1.
sigmas(0) = 0.
Gamma0 = cj * omega * SQRT(murs(0)*epsilonrs(0))/c0
Gammas(0) = Gamma0
thetas(0) = CMPLX(theta)
statc = 0
psi(0) = (0.,0.)
DO i = 1, n_lay
    murs(i) = Materials( INT(WTypes(wt,(1+n_lay+i))) , 2) ! check it
    epsilonrs(i) = Materials( INT(WTypes(wt,(1+n_lay+i))) , 1)
    sigmas(i) = Materials( INT(WTypes(wt,(1+n_lay+i))) , 3) / 1000.
    Gammas(i) = SQRT( cj*omega*murs(i)*mu0 * (sigmas(i) + cj*omega*epsilonrs(i)*epsilon0) )

    temp = ( Gammas(i-1)/Gammas(i)*CSIN(thetas(i-1)) )
    
    IF (temp == (0.,0.)) THEN
       thetas(i) = thetas(i-1)
    ELSE
    temp = 1./temp
    !thetas(i) = -cj*CLOG(cj*temp+CSQRT(temp**2-1.)) 
     thetas(i) = -cj*CLOG( CSQRT(1.-1./temp**2) + cj/temp )
    END IF
    
    psi(i) = WTypes(wt, 1+i)/100.*Gammas(i)*CCOS(thetas(i)) ! psi : 1:n_lay
    loss_tans(i) = sigmas(i)/(omega*epsilon0*epsilonrs(i))
END DO
Gammas(n_lay+1) = CMPLX(Gamma0)
temp = ( Gammas(n_lay)/Gammas(n_lay+1)*CSIN(thetas(n_lay)) )
! loss_tans, murs, epsilonrs, sigmas
murs(n_lay+1) = mu0
epsilonrs(n_lay+1) = 1.
sigmas(n_lay+1) = 0.
loss_tans(n_lay+1) = 0.
IF (temp == (0.,0.)) THEN
       thetas(n_lay+1) = thetas(n_lay)
    ELSE
    temp = 1./temp
    !thetas(i) = -cj*CLOG(cj*temp+CSQRT(temp**2-1.)) 
     thetas(n_lay+1) = -cj*CLOG( CSQRT(1.-1./temp**2) + cj/temp )
    END IF

!thetas(n_lay+1) = -cj*CLOG(cj*temp+CSQRT(1.-temp**2)) 
!psi(n_lay+1) = WTypes(wt, 1+n_lay+1)/100*Gammas(n_lay+1)*CCOS(thetas(n_lay+1))
A(n_lay+1) = CMPLX(1.) ! check these all
D(n_lay+1) = CMPLX(1.)
B(n_lay+1) = CMPLX(0.)
 C(n_lay+1) = CMPLX(0.)
loss_tans(n_lay+1) = 0.
Y(n_lay+1) = CCOS(thetas(n_lay+1))/CCOS(thetas(n_lay)) * &
SQRT( epsilon0*epsilonrs(n_lay+1)*(1.-cj*loss_tans(n_lay+1)) / epsilon0/epsilonrs(n_lay)/(1.-cj*loss_tans(n_lay)) )
Z(n_lay+1) = CCOS(thetas(n_lay+1))/CCOS(thetas(n_lay)) * &
SQRT( epsilon0*epsilonrs(n_lay)*(1.-cj*loss_tans(n_lay)) / epsilon0/epsilonrs(n_lay+1)/(1.-cj*loss_tans(n_lay+1)) )

ptemp = (0.,0.)
mtemp = (0.,0.)
layersloop: DO j=n_lay,0,-1
    IF (j>=0) THEN
       Y(j+1) = CCOS(thetas(j+1))/CCOS(thetas(j+1-1)) * &
       SQRT( epsilon0*epsilonrs(j+1)*(1.-cj*loss_tans(j+1)) / epsilon0/epsilonrs(j+1-1)/(1.-cj*loss_tans(j+1-1)) )
       Z(j+1) = CCOS(thetas(j+1))/CCOS(thetas(j+1-1)) * & 
       SQRT( epsilon0*epsilonrs(j+1-1)*(1.-cj*loss_tans(j+1-1)) / epsilon0/epsilonrs(j+1)/(1.-cj*loss_tans(j+1)) )
    END IF
    
    !ptemp = ptemp + psi(j)
    !mtemp = mtemp - psi(j)
   ! WRITE(*,*) 'psi(',j,') = ', psi(j)
    !A(j) = 1./2.*(A(j+1)*(1.+Y(j+1))+B(j+1)*(1-Y(j+1)))
    !B(j) = 1./2.*(A(j+1)*(1.-Y(j+1))+B(j+1)*(1+Y(j+1)))
    !C(j) = 1./2.*  (  C(j+1)*(1.+Z(j+1))  +  D(j+1)*(1-Z(j+1))  )
    !D(j) = 1./2.*(C(j+1)*(1.-Z(j+1))+D(j+1)*(1+Z(j+1)))
    IF (ABS(REAL(psi(j)))>=87.) THEN
            statc = 1
    END IF
    
    A(j) = CEXP(psi(j))/2.*(A(j+1)*(1.+Y(j+1))+B(j+1)*(1-Y(j+1)))
    B(j) = CEXP(-psi(j))/2.*(A(j+1)*(1.-Y(j+1))+B(j+1)*(1+Y(j+1)))
    C(j) = CEXP(-psi(j))/2.*  (  C(j+1)*(1.+Z(j+1))  +  D(j+1)*(1-Z(j+1))  )
    D(j) = CEXP(psi(j))/2.*(C(j+1)*(1.-Z(j+1))+D(j+1)*(1+Z(j+1)))
    
    IF ( ((A(j)-A(j)) /= CMPLX(0.)) .AND. ((D(j)-D(j))/=CMPLX(0.)) ) THEN
           statc = 1
          EXIT layersloop
  END IF 
END DO layersloop
IF ((statc == 0) .AND. (D(0)/=CMPLX(0.)) .AND. (A(0)/=CMPLX(0.)) ) THEN
Rpar = C(0)/D(0)
Rnorm = B(0)/A(0)
!IF (Rpar/=Rpar) THEN
!        WRITE(*,*) 'Rpar coef error'
!        WRITE(*,*) 'Rpar coef error'
!END IF
!IF (Rnorm/=Rnorm) THEN
!        WRITE(*,*) 'Rnorm coef error'
!        WRITE(*,*) 'Rnorm coef error'
!END IF
Tpar = 1./D(0)
Tnorm = 1./A(0)
ELSE
        STOP 'ERROR: material too wide and too conductive'
END IF


RETURN
END SUBROUTINE R_T_coeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION E_first_hit(loc_hit,P,f)
!USE ALL_DATA
IMPLICIT NONE
REAL, DIMENSION(3),INTENT(IN) :: loc_hit
COMPLEX, DIMENSION(3):: E_first_hit
REAL :: D, r_dist, theta, phi, P,f
COMPLEX :: EE

P = TX(4)
 CALL Get_angle(loc_hit, r_dist, theta, phi)
D = Directivity(theta, phi, 1)
EE = cj * SQRT(P/36.5/1000.) / (2.*pi*epsilon0*c0)*D

E_first_hit = (/ EE*COS(theta)*COS(phi) , EE*COS(theta)*SIN(phi), -EE*SIN(theta)  /)

RETURN
END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Path_length(Paths_sbrtn, Path_lengths)
!USE ALL_DATA
IMPLICIT NONE
REAL, DIMENSION(4,R+3),INTENT(IN) :: Paths_sbrtn
REAL, INTENT(OUT) :: Path_lengths
REAL :: temp
INTEGER i, j
REAL, DIMENSION(3) :: first_point, second_point
! NOT very good ... null rows of All_Paths_sbrtn should be deleted.

     first_point = Paths_sbrtn ( 2:4,2)
     temp = 0.
     j = 2
     each_path : DO
                    j = j+1
                    second_point = Paths_sbrtn( 2:4,j)
                    temp = temp + SQRT( (first_point(1)-second_point(1))**2&
                    +(first_point(2)-second_point(2))**2+(first_point(3)-second_point(3))**2 )
                    
                    IF (Paths_sbrtn(1,j)==0.) EXIT each_path
                    first_point = second_point
                 END DO each_path
     
     Path_lengths = temp
!Good_paths = All_Paths_sbrtn(1:num, :, :)
RETURN
END SUBROUTINE Path_length
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION Path_field( path , pl )
         !USE ALL_DATA
          IMPLICIT NONE
         REAL, DIMENSION(4,R+3), INTENT(IN):: path
         REAL, INTENT(IN) :: pl
         COMPLEX, DIMENSION(3) :: Path_field
          INTEGER, DIMENSION(2) :: sz
          INTEGER :: lp, ORD, i, w_loc_num, w_num, loctr3,ii
          REAL, DIMENSION(:,:), ALLOCATABLE :: thispath
          REAL, DIMENSION(3) :: loc_obs, loc_hit, ploc_tx
          COMPLEX, DIMENSION(3) :: Ei
          REAL :: r_dist, theta, phi, D, theta_inc
          REAL :: Beta, freq
          COMPLEX :: Rpar, Rnorm, Tpar, Tnorm
          Beta = 2.*pi*TX(5)*SQRT(epsilon0*mu0)
          freq = TX(5)
          
         ! Other declarations
         sz = SHAPE (path)
         lp = sz (2)
         ORD = INT(path(1,1))
         ! path : 4 X (R+3)
         Allocate( thispath (4, (ORD+3) ) )
          thispath = path(:,1:ORD+3)
         loc_obs = thispath(2:4,ORD+2) !gets the first incident 
          !which might be obs if ORD=0.
         Ei = E_first_hit(loc_obs, TX(4) , TX(5))
         !CALL Get_angle(loc_obs, r_dist, theta, phi)
          
         !D = Directivity(theta,phi,1)
         !Ei = D * Ei
         ploc_tx = loc_tx
       IF (ORD>0) THEN  
       DO i = 1, ORD
        w_loc_num = thispath(1,ORD+3-i)
        w_num = INT(Tree(0, w_loc_num))
        loc_hit = thispath( 2:4 , ORD+3-i)
        !IF (loc_hit(1) == 4. .AND. loc_hit(2)==5. AND. ploc_tx(1)==4. AND. ploc_tx(2)==5. ) THEN
        !IF ( loc_hit(2)==5.  .AND. ploc_tx(2)==5. ) THEN
        !        DO ii=1,4
        !                WRITE(*,*) thispath(ii,:)
        !        END DO
!
 !       END IF


        CALL get_inc_theta(ploc_tx, loc_hit, w_num, theta_inc)
        !WRITE(*,*)'theta_inc = ',theta_inc
        !WRITE(*,*) 'the INT thing = ', INT(WTypes(All_Walls(w_num,15),1))
        !WRITE(*,*) ' w_num = ', w_num
!        WRITE(*,*) 'theta_inc(',i,')= ',theta_inc
        CALL R_T_coeffs( theta_inc, freq, w_num, &
        INT(WTypes(INT(Walls(w_num,7)),1)), Rpar, Rnorm, Tpar, Tnorm)
        !WRITE(*,*) 'Rpar(',i,')= ', (Rpar)
        !WRITE(*,*) 'Rnorm(',i,')= ',(Rnorm)
        !IF (Ei(1)/=Ei(1) .OR. Ei(2)/=Ei(2) .OR. Ei(3)/=Ei(3)) &
        !STOP 'Error before calling Rebuild_E'
        CALL Rebuild_E( Ei, w_loc_num, Rpar, Rnorm, ploc_tx, loc_hit)
        !IF (Ei(1)/=Ei(1) .OR. Ei(2)/=Ei(2) .OR. Ei(3)/=Ei(3)) &
        !STOP 'Error after calling Rebuild_E'
 !       WRITE(*,*)'Ei after(',i,')= ',Ei
        ploc_tx = loc_hit
      END DO
      END IF
      IF ((Beta*pl) > 8E5) WRITE(*,*) 'pl = ', pl
      Path_field = Ei * CEXP(-cj * Beta * pl) / pl
      RETURN
      END FUNCTION Path_field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_inc_theta(ploc_tx, loc_hit, w_num3, theta_inc)
!USE ALL_DATA
IMPLICIT NONE
REAL, DIMENSION(3), INTENT(IN)::ploc_tx
REAL, DIMENSION(3), INTENT(IN)::loc_hit
INTEGER, INTENT(IN) :: w_num3
REAL, INTENT(OUT) :: theta_inc
REAL :: n_vec, mytempval
REAL, DIMENSION(3) :: minus_hit , vec_n_vec
INTEGER :: theta_stat
!WRITE(*,*) 'ploc = ',ploc_tx
!WRITE(*,*) 'loc_hit = ',loc_hit
theta_stat = 0
IF (All_Walls(w_num3,13) == 1.) THEN
        !mytempval = ploc_tx(1)-loc_hit(1)
        n_vec = (ploc_tx(1)-loc_hit(1))/ABS((ploc_tx(1)-loc_hit(1)))
        IF (n_vec /= n_vec) theta_stat = 1
        !IF (mytempval == 0.) THEN
               ! vec_n_vec =(
        vec_n_vec = (/ n_vec, 0. , 0. /)
        !write(*,*) 'hello 1'
ELSE IF (All_Walls(w_num3,13) == 2.) THEN
       n_vec = (ploc_tx(2)-loc_hit(2))/ABS((ploc_tx(2)-loc_hit(2)))
        IF (n_vec /= n_vec) theta_stat = 1
       vec_n_vec = (/ 0.,n_vec, 0. /)
        !write(*,*) 'hello 2'
        !write(*,*) 'vals =', ploc_tx(2),loc_hit(2)
        !write(*,*) 'n_vec= ', n_vec
ELSE IF (All_Walls(w_num3,13) == 3.) THEN
       n_vec = (ploc_tx(3)-loc_hit(3))/ABS((ploc_tx(3)-loc_hit(3)))
        IF (n_vec /= n_vec) theta_stat = 1
       vec_n_vec = (/ 0., 0. , n_vec /)
       !write(*,*) 'hello 3'
END IF
minus_hit = ploc_tx - loc_hit
!WRITE(*,*) 'hello',(DOT_PRODUCT(vec_n_vec,vec_n_vec)),(DOT_PRODUCT(minus_hit,minus_hit)) 
!IF (( ACOS( DOT_PRODUCT(vec_n_vec,minus_hit)/ &
!   SQRT( (DOT_PRODUCT(vec_n_vec,vec_n_vec))*(DOT_PRODUCT(minus_hit,minus_hit)) ) )) /= & 
!   ( ACOS( DOT_PRODUCT(vec_n_vec,minus_hit)/ &
!   SQRT( (DOT_PRODUCT(vec_n_vec,vec_n_vec))*(DOT_PRODUCT(minus_hit,minus_hit)) ) ))  ) THEN
        !WRITE(*,*) 'DOT_PRODUCT(vec_n_vec,minus_hit)=', DOT_PRODUCT(vec_n_vec,minus_hit) 
        !WRITE(*,*) 'SQRT( (DOT_PRODUCT(vec_n_vec,vec_n_vec))*(DOT_PRODUCT(minus_hit,minus_hit)) ) = ',SQRT( (DOT_PRODUCT(vec_n_vec,vec_n_vec))*(DOT_PRODUCT(minus_hit,minus_hit)) ) 
!END IF
!WRITE(*,*) 'MINUS_HIT = ', minus_hit
!WRITE(*,*) 'ploc_tx) = ',ploc_tx
!WRITE(*,*) 'loc_hit = ', loc_hit
IF (theta_stat == 0) THEN
theta_inc = ACOS (DOT_PRODUCT(vec_n_vec,minus_hit)/ &
SQRT( (DOT_PRODUCT(vec_n_vec,vec_n_vec))*(DOT_PRODUCT(minus_hit,minus_hit)) ) )
ELSE
        theta_inc = 0.
END IF
!WRITE(*,*) 'theta_inc = ', theta_inc*180./pi
RETURN
END SUBROUTINE get_inc_theta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Rebuild_E( Ei, w_loc_num, Rpar, Rnorm, ploc_tx, loc_hit)
!USE ALL_DATA
IMPLICIT NONE
 COMPLEX, DIMENSION(3), INTENT(INOUT) :: Ei
 INTEGER, INTENT(IN) :: w_loc_num
 COMPLEX, INTENT(IN) :: Rpar, Rnorm
 REAL, DIMENSION(3), INTENT(IN) :: ploc_tx, loc_hit

 REAL, DIMENSION(3) :: norm2pln_inc,a,b,n, minus_hit, vec_n_vec
 COMPLEX, DIMENSION(3) :: norm_part, par_part, par_part2, Eii
 REAL :: d, n_vec
 INTEGER i, w_num5, walln, statinc

 w_num5 = INT(Tree(0,w_loc_num))
 !WRITE(*,*) 'd = ',d
statinc = 0
IF (All_Walls(w_num5,13) == 1.) THEN
       n_vec = (ploc_tx(1)-loc_hit(1))/ABS((ploc_tx(1)-loc_hit(1)))
       vec_n_vec = (/ n_vec, 0. , 0. /)
       walln = 1
ELSE IF (All_Walls(w_num5,13) == 2.) THEN
       n_vec = (ploc_tx(2)-loc_hit(2))/ABS((ploc_tx(2)-loc_hit(2)))
       vec_n_vec = (/ 0.,n_vec, 0. /)
       walln = 2
ELSE IF (All_Walls(w_num5,13) == 3.) THEN
       n_vec = (ploc_tx(3)-loc_hit(3))/ABS((ploc_tx(3)-loc_hit(3)))
       vec_n_vec = (/ 0., 0. , n_vec /)
       walln = 3
ELSE
       READ(*,*)
END IF
! IF (ABS(d)==1.) THEN
!       n = (/ d/ABS(d) , 0. , 0. /)
! ELSE IF (ABS(d)==2.) THEN
!       n = (/ 0., d/ABS(d) , 0. /)
! ELSE IF (ABS(d)==3.) THEN
!       n = (/ 0., 0. , d/ABS(d) /)
! END IF
!WRITE(*,*) 'w_loc_num = ', w_loc_num
!WRITE(*,*) 'w_num) = ', w_num5
!WRITE(*,*) vec_n_vec
 n = vec_n_vec
 minus_hit = -loc_hit + ploc_tx
 a = n
 b = minus_hit
 !WRITE(*,*) 'a = ',a
 !WRITE(*,*) 'b = ',b
 norm2pln_inc = (/ a(2)*b(3)-a(3)*b(2), -a(1)*b(3)+a(3)*b(1), a(1)*b(2)-a(2)*b(1) /)
 !WRITE(*,*) 'norm2pln_inc = ',norm2pln_inc
 IF ((DOT_PRODUCT(norm2pln_inc,norm2pln_inc))==0. .OR. &
 norm2pln_inc(1)/=norm2pln_inc(1) .OR. norm2pln_inc(2)/=norm2pln_inc(2)&
  .OR. norm2pln_inc(3)/=norm2pln_inc(3) )THEN
         statinc = 1
  !b= b-(/0.1,0.1,0.1/)
  !norm2pln_inc = (/ a(2)*b(3)-a(3)*b(2), -a(1)*b(3)+a(3)*b(1), a(1)*b(2)-a(2)*b(1) /)
END IF
IF (statinc == 0) THEN
 norm2pln_inc = norm2pln_inc / SQRT(DOT_PRODUCT(norm2pln_inc,norm2pln_inc))
Eii = Ei
 norm_part = CMPLX(1.) * norm2pln_inc * &
 (  Eii(1)*norm2pln_inc(1) + Eii(2)*norm2pln_inc(2) + Eii(3)*norm2pln_inc(3)  )
 par_part = Eii - norm_part
 par_part2 = par_part
 par_part2(walln) =  par_part(walln) * (-1.)
 par_part2 = -par_part2

 Ei = norm_part * Rnorm + par_part2 * Rpar
ELSE 
        Ei = Ei*Rnorm
END IF
 RETURN
 END SUBROUTINE Rebuild_E
!888888888888888888888888888888888888
!888888888888888888888888888888888888
!888888888888888888888888888888888888888888888
!888888888888888888888888888888888888888888888
!888888888888888888888888888888888888888888888888
!888888888888888888888888888888888888888888888
!88888888888888888888888888888888888888888888
!88888888888888888888888888888888888888888888
 
 END MODULE ALL_DATA


!888888888888888888888888888888888888888888888888888

!88888888888888888888888888888888888888888888888888888
!88888888888888888888888888888888888888888888888888888
!88888888888888888888888888888888888888888888888888888

! code starts
MODULE commandresults
!!!!!!!!!  
  USE ALL_DATA

!!!!!!!!!!
 CONTAINS
!888888888888888888888888888888888
SUBROUTINE GET_WORD(the_word,nn)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nn
 CHARACTER(len=nn), INTENT(OUT) :: the_word
 INTEGER i
 
 ! ignoring blanks
 delete_blanks1: DO
        IF (line1(1:1) /= ' ') EXIT delete_blanks1
        line1 = line1(2:)
 END DO delete_blanks1
 ! setting "the_word" to blanks
 DO i=1,nn
        the_word(i:i)=' '
 END DO
 i=1
 ! get "the_word". stop at blank
 word: DO
       the_word(i:i)=line1(1:1)
        line1 = line1(2:)
        i=i+1
        IF (line1(1:1) == ' ') EXIT word
        
 END DO word
 !WRITE(*,*)'from GET_WORD, the_word is :',the_word
 !
 !the_word = line1(1:nn)
 !
  ! ignoring blanks again
 delete_blanks2: DO
        IF (line1(1:1) /= ' ') EXIT delete_blanks2
        line1 = line1(2:)
 END DO delete_blanks2
 RETURN
END SUBROUTINE GET_WORD
!888888888888888888888888888888888
SUBROUTINE GET_WORD_var(the_word,nn)
 IMPLICIT NONE
 INTEGER, INTENT(OUT) :: nn
 CHARACTER(len=10), INTENT(OUT) :: the_word
 INTEGER i,j
 ! ignoring blanks
 delete_blanks1: DO
        IF (line1(1:1) /= ' ') EXIT delete_blanks1
        line1 = line1(2:)
 END DO delete_blanks1
 ! setting "the_word" to blanks
 DO i=1,10
        the_word(i:i)=' '
 END DO
 j=1
 ! get "the_word". stop at blank
 word: DO
       the_word(j:j)=line1(1:1)
        j=j+1
        line1 = line1(2:)
        IF (line1(1:1) == ' ') EXIT word
        
 END DO word
 nn=j-1
 !the_word = line1(1:nn)
 !
  ! ignoring blanks again
 j = 0
 delete_blanks2: DO
       j= j+1
        !write(*,*) 'inside loop from inside GET_WORD_var, line1 :',line1,'x'
       !write(*,*) 'LEN(line1) =:',len(line1)
       IF (j>LEN(line1)) EXIT delete_blanks2
        IF (line1(1:1) /= ' ') EXIT delete_blanks2
        line1 = line1(2:)
 END DO delete_blanks2
 RETURN
END SUBROUTINE GET_WORD_var
!888888888888888888888888888888888
SUBROUTINE COM_MATERIAL(nums)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nums  
  REAL :: epsr,sigma,mur
  CHARACTER(len=10)::mat_name

  CALL GET_WORD(mat_name,10)
  material_names(nums) = mat_name
  READ(line1,*) epsr,sigma,mur
  !WRITE(*,*) 'material name :',material_names(nums)
  !WRITE(*,*) 'epsr =',epsr,'  sigma =',sigma,'   mur =',mur
  material_mat(nums,:) = (/ epsr, sigma,mur/)
END SUBROUTINE COM_MATERIAL
!888888888888888888888888888888888
SUBROUTINE COM_WALL_TYPE(nums)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nums
  REAL :: a, wll_wdth,g
  INTEGER:: name_len, i, j,laynum
  CHARACTER(len=10)::wlltp_name, temp, wlltp_sl,my_mat_name
  ! get wall_type_name: blah blah
  CALL GET_WORD_var(wlltp_name,name_len)
  wall_type_names(nums) = wlltp_name
  ! get SOLID or LAYERED

  CALL GET_WORD_var(wlltp_sl,name_len) !sl:SOLID or LAYERED
  IF (name_len == 5 .AND. wlltp_sl == 'SOLID') THEN
         ! get wall material's name
        CALL GET_WORD_var(my_mat_name,name_len)
       ! find which name in "material_names" matches "my_mat_name"
        DO i=1,8 ! 8 is max. number of materials
              temp = material_names(i)
              IF (temp(1:name_len)==my_mat_name) g = REAL(i)
        END DO
        READ(line1,*) wll_wdth
         wall_type_mat(nums,:) = (/1. ,g, wll_wdth ,0.,0. ,0. , 0. ,0. ,0. ,0. ,0. ,0. ,0. ,0. ,0./)
  ELSEIF (name_len == 7 .AND. wlltp_sl == 'LAYERED') THEN
       ! get the laynum number
       laynum = IACHAR(line1(1:1))-IACHAR('0')
       wall_type_mat(nums,1)=laynum
       ! get rid of that character and blanks
       line1 = line1(2:)
        ! get material name and width laynum times
       DO j=1,laynum
              CALL GET_WORD_var(my_mat_name,name_len)
              ! find which name in "material_names" matches "my_mat_name"
               DO i=1,8 ! 8 is max. number of materials
                     temp = material_names(i)
                     IF (temp(1:name_len)==my_mat_name) g = REAL(i)
               END DO
              ! getting the width in char mode
               CALL GET_WORD_var(my_mat_name,name_len)!my_mat_name is inappropriate but never mind!
              ! converting width from char to read
              READ(my_mat_name,*) wll_wdth
              wall_type_mat(nums,(1+(j-1)*2+1):(1+(j-1)*2+2)) = (/g, wll_wdth /)
              
       END DO
       !write(*,*)'nums,walls :',nums,wall_type_mat(nums,:)

  END IF
  
 RETURN
END SUBROUTINE COM_WALL_TYPE
!888888888888888888888888888888888
SUBROUTINE COM_POINT(nums)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nums
  INTEGER :: i,j,ln
  CHARACTER(len=10):: tempchar,pnt_nm
  REAL :: x,y,z
  ! setting point_names(nums) to all blanks
  DO i=1,10
       tempchar(i:i)=' '
  ENDDO
  point_names(nums) = tempchar
  !getting the name of the point
  CALL GET_WORD_var(pnt_nm,ln)
  point_names(nums) = pnt_nm
  READ(line1,*) x,y,z
  point_mat(nums,:) = (/x,y,z/)
  !WRITE(*,*) 'point name :',point_names(nums)
  !WRITE(*,*) 'point coordinates :',point_mat(nums,:)

END SUBROUTINE COM_POINT
!888888888888888888888888888888888
SUBROUTINE COM_DEFINE_WALL(w_nums,p_nums)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: w_nums, p_nums
  INTEGER :: ln, i,j
  CHARACTER(len=10)::pnt_nm,tempchar
  tempchar = '          '
  !get the first point name
  CALL GET_WORD_var(pnt_nm,ln)
  ! find the first point name in stack
  findpntname1: DO  i=1,p_nums
         IF (point_names(i) == pnt_nm) THEN
                define_wall_mat(w_nums,1)=i
                EXIT findpntname1
         ENDIF
  ENDDO findpntname1
  !get the second point name
  CALL GET_WORD_var(pnt_nm,ln)
  ! find the second point name in stack
  findpntname2: DO  i=1,p_nums
         IF (point_names(i) == pnt_nm) THEN
                define_wall_mat(w_nums,2)=i
                EXIT findpntname2
         ENDIF
  ENDDO findpntname2
  ! get the wall_type name
  CALL GET_WORD_var(tempchar,ln)
  ! find the the wall_type in stack
  findwall: DO i=1,8 ! 8 is the max. no. of wall types
         IF (wall_type_names(i)==tempchar) THEN
              define_wall_types(w_nums) = i
              EXIT findwall
       ENDIF
  ENDDO findwall
  !report
  !WRITE(*,*)':-D ','DFEINE_WALL ',point_names(define_wall_mat(w_nums,1)),'  ',point_names(define_wall_mat(w_nums,2)),'  ',wall_type_names(define_wall_types(w_nums))
END SUBROUTINE COM_DEFINE_WALL
!888888888888888888888888888888888
SUBROUTINE COM_FREQUENCY()
  IMPLICIT NONE
  REAL:: tempfreq
  INTEGER :: ln, i,j
  CHARACTER(len=10)::pnt_nm,tempchar
  CALL GET_WORD_var(tempchar,ln)
  !WRITE(*,*) 'from COM_FREQUENCY, line :',line1
  READ(tempchar,*) tempfreq
  TX(5)=tempfreq*1.E6
  CALL GET_WORD_var(tempchar,ln)
  IF (.NOT.(tempchar=='MHz' .OR. tempchar=='MHz.')) THEN
         WRITE(*,*)'only "MHz"'
       STOP
  ENDIF
END SUBROUTINE COM_FREQUENCY
!888888888888888888888888888888888
SUBROUTINE COM_SOURCE_LOCATION(pnt_num)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: pnt_num
  INTEGER :: ln, i,jtemp
  CHARACTER(len=10)::tempchar
  CALL GET_WORD_var(tempchar,ln)
  !finding the point
  fnd_pnt: DO i=1,pnt_num
     IF (point_names(i)==tempchar) THEN
       jtemp=i
       EXIT fnd_pnt
     ENDIF
  ENDDO fnd_pnt
  loc_tx = point_mat(jtemp,:)
  TX(1:3) = loc_tx
END SUBROUTINE COM_SOURCE_LOCATION
!888888888888888888888888888888888
SUBROUTINE COM_SOURCE_POWER()
  IMPLICIT NONE
  INTEGER :: ln, i,jtemp
  CHARACTER(len=10)::tempchar
  REAL :: tempower
  CALL GET_WORD_var(tempchar,ln)
  READ(tempchar,*) tempower
  TX(4) = tempower
  CALL GET_WORD_var(tempchar,ln)
  IF (.NOT.(tempchar=='mW' .OR. tempchar=='mW.')) THEN
         WRITE(*,*)'only "mW"'
       STOP
  ENDIF
END SUBROUTINE COM_SOURCE_POWER
!888888888888888888888888888888888
SUBROUTINE COM_DIPOLE()
  IMPLICIT NONE
  INTEGER :: ln, i,jtemp

  REAL :: x,y,z

  READ(line1,*) x,y,z
  WRITE(*,*)'OK since only DIPOLE is accepted.'
  IF (.NOT.(x==0. .AND. y==0. .AND. z==1.)) THEN
         WRITE(*,*)'DIPOLE only in z direction'
       STOP
  ENDIF
  TX(6)=1.
END SUBROUTINE COM_DIPOLE
!888888888888888888888888888888888
SUBROUTINE COM_RCVRPOINT(oneline)
  IMPLICIT NONE
  CHARACTER(len=80), INTENT(IN) :: oneline
END SUBROUTINE COM_RCVRPOINT
!888888888888888888888888888888888
SUBROUTINE COM_RCVRLINE(pnt_nm)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: pnt_nm
  INTEGER :: ln, i,jtemp
  CHARACTER(len=10)::tempchar
  REAL, DIMENSION(3):: start_pnt,end_pnt
  CALL GET_WORD_var(tempchar,ln)
   !finding the start point
  fnd_pnt: DO i=1,pnt_nm
     IF (point_names(i)==tempchar) THEN
       jtemp=i
       EXIT fnd_pnt
     ENDIF
  ENDDO fnd_pnt
  start_pnt = point_mat(jtemp,:)
  CALL GET_WORD_var(tempchar,ln)
   !finding the end point
  fnd_pnt2: DO i=1,pnt_nm
     IF (point_names(i)==tempchar) THEN
       jtemp=i
       EXIT fnd_pnt2
     ENDIF
  ENDDO fnd_pnt2
  end_pnt = point_mat(jtemp,:)
  !finding xy(3)?,xz(2)?,yz(1)?
  jtemp = 0
  DO i=1,3
         IF (start_pnt(i)==end_pnt(i)) THEN
              area_indicator=1 !This is now actually line indicator. along x,y or
              !axes only
              jtemp = jtemp + 1
       ENDIF
  ENDDO
  IF (jtemp /= 2) THEN
       WRITE(*,*) 'Next time line must be along either x, y or z axis. Good boy!'
       ! STOP
  ENDIF
  READ(line1,*) delta_incr
  !WRITE(*,*) 'start =',start_pnt
  !WRITE(*,*) 'end   =',end_pnt
  !WRITE(*,*) 'indicator =',area_indicator
  !WRITE(*,*) 'delta_incr =', delta_incr
  xn = INT(ABS(end_pnt(1)-start_pnt(1))/(delta_incr/100.))
  yn = INT(ABS(end_pnt(2)-start_pnt(2))/(delta_incr/100.))
  zn = INT(ABS(end_pnt(3)-start_pnt(3))/(delta_incr/100.))
  ALLOCATE(xrec(xn+1))
  ALLOCATE(yrec(yn+1))
  ALLOCATE(zrec(zn+1))
  DO i=1,xn+1
       xrec(i)=MINVAL((/end_pnt(1),start_pnt(1) /))+(i-1)*delta_incr/100.
  ENDDO
  DO i=1,yn+1
       yrec(i)=MINVAL((/end_pnt(2),start_pnt(2) /))+(i-1)*delta_incr/100.
  ENDDO
  DO i=1,zn+1
       zrec(i)=MINVAL((/end_pnt(3),start_pnt(3) /))+(i-1)*delta_incr/100.
  ENDDO
  !IF (area_indicator==1) THEN
  !     ALLOCATE(Ex_abs_rec(yn+1,zn+1))
  !     ALLOCATE(Ex_phase_rec(yn+1,zn+1))
  !     ALLOCATE(Ey_abs_rec(yn+1,zn+1))
  !     ALLOCATE(Ey_phase_rec(yn+1,zn+1))
  !     ALLOCATE(Ez_abs_rec(yn+1,zn+1))
  !     ALLOCATE(Ez_phase_rec(yn+1,zn+1))
  !     ALLOCATE(Et_abs_rec(yn+1,zn+1))
  !     ALLOCATE(Emean_rec(yn+1,zn+1))
  !ELSEIF (area_indicator==2) THEN
  !     ALLOCATE(Ex_abs_rec(xn+1,zn+1))
  !     ALLOCATE(Ex_phase_rec(xn+1,zn+1))
  !     ALLOCATE(Ey_abs_rec(xn+1,zn+1))
  !     ALLOCATE(Ey_phase_rec(xn+1,zn+1))
  !     ALLOCATE(Ez_abs_rec(xn+1,zn+1))
  !     ALLOCATE(Ez_phase_rec(xn+1,zn+1))
  !     ALLOCATE(Et_abs_rec(xn+1,zn+1))
  !     ALLOCATE(Emean_rec(xn+1,zn+1))
  !ELSEIF (area_indicator==3) THEN
  !     ALLOCATE(Ex_abs_rec(xn+1,yn+1))
  !     ALLOCATE(Ex_phase_rec(xn+1,yn+1))
  !     ALLOCATE(Ey_abs_rec(xn+1,yn+1))
  !     ALLOCATE(Ey_phase_rec(xn+1,yn+1))
  !     ALLOCATE(Ez_abs_rec(xn+1,yn+1))
  !     ALLOCATE(Ez_phase_rec(xn+1,yn+1))
  !     ALLOCATE(Et_abs_rec(xn+1,yn+1))
  !     ALLOCATE(Emean_rec(xn+1,yn+1))
  !ENDIF
END SUBROUTINE COM_RCVRLINE
!888888888888888888888888888888888
SUBROUTINE COM_RCVRGRID(pnt_nm)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: pnt_nm
  INTEGER :: ln, i,jtemp
  CHARACTER(len=10)::tempchar
  REAL, DIMENSION(3):: start_pnt,end_pnt
  CALL GET_WORD_var(tempchar,ln)
   !finding the start point
  fnd_pnt: DO i=1,pnt_nm
     IF (point_names(i)==tempchar) THEN
       jtemp=i
       EXIT fnd_pnt
     ENDIF
  ENDDO fnd_pnt
  start_pnt = point_mat(jtemp,:)
  CALL GET_WORD_var(tempchar,ln)
   !finding the end point
  fnd_pnt2: DO i=1,pnt_nm
     IF (point_names(i)==tempchar) THEN
       jtemp=i
       EXIT fnd_pnt2
     ENDIF
  ENDDO fnd_pnt2
  end_pnt = point_mat(jtemp,:)
  !finding xy(3)?,xz(2)?,yz(1)?
  DO i=1,3
         IF (start_pnt(i)==end_pnt(i)) THEN
              area_indicator=i
              jtemp=0
       ENDIF
  ENDDO
  IF (jtemp /= 0) THEN
       WRITE(*,*) 'area on either xy,xz,yz planes only'
        STOP
  ENDIF
  READ(line1,*) delta_incr
  !WRITE(*,*) 'start =',start_pnt
  !WRITE(*,*) 'end   =',end_pnt
  !WRITE(*,*) 'indicator =',area_indicator
  !WRITE(*,*) 'delta_incr =', delta_incr
  xn = INT(ABS(end_pnt(1)-start_pnt(1))/(delta_incr/100.))
  yn = INT(ABS(end_pnt(2)-start_pnt(2))/(delta_incr/100.))
  zn = INT(ABS(end_pnt(3)-start_pnt(3))/(delta_incr/100.))
  ALLOCATE(xrec(xn+1))
  ALLOCATE(yrec(yn+1))
  ALLOCATE(zrec(zn+1))
  DO i=1,xn+1
       xrec(i)=MINVAL((/end_pnt(1),start_pnt(1) /))+(i-1)*delta_incr/100.
  ENDDO
  DO i=1,yn+1
       yrec(i)=MINVAL((/end_pnt(2),start_pnt(2) /))+(i-1)*delta_incr/100.
  ENDDO
  DO i=1,zn+1
       zrec(i)=MINVAL((/end_pnt(3),start_pnt(3) /))+(i-1)*delta_incr/100.
  ENDDO
 ! IF (area_indicator==1) THEN
 !      ALLOCATE(Ex_abs_rec(yn+1,zn+1))
 !      ALLOCATE(Ex_phase_rec(yn+1,zn+1))
 !      ALLOCATE(Ey_abs_rec(yn+1,zn+1))
 !      ALLOCATE(Ey_phase_rec(yn+1,zn+1))
 !      ALLOCATE(Ez_abs_rec(yn+1,zn+1))
 !      ALLOCATE(Ez_phase_rec(yn+1,zn+1))
 !      ALLOCATE(Et_abs_rec(yn+1,zn+1))
 !      ALLOCATE(Emean_rec(yn+1,zn+1))
 ! ELSEIF (area_indicator==2) THEN
 !      ALLOCATE(Ex_abs_rec(xn+1,zn+1))
 !      ALLOCATE(Ex_phase_rec(xn+1,zn+1))
 !      ALLOCATE(Ey_abs_rec(xn+1,zn+1))
 !      ALLOCATE(Ey_phase_rec(xn+1,zn+1))
 !      ALLOCATE(Ez_abs_rec(xn+1,zn+1))
 !      ALLOCATE(Ez_phase_rec(xn+1,zn+1))
 !      ALLOCATE(Et_abs_rec(xn+1,zn+1))
 !      ALLOCATE(Emean_rec(xn+1,zn+1))
 ! ELSEIF (area_indicator==3) THEN
 !      ALLOCATE(Ex_abs_rec(xn+1,yn+1))
 !      ALLOCATE(Ex_phase_rec(xn+1,yn+1))
 !      ALLOCATE(Ey_abs_rec(xn+1,yn+1))
 !      ALLOCATE(Ey_phase_rec(xn+1,yn+1))
 !      ALLOCATE(Ez_abs_rec(xn+1,yn+1))
 !      ALLOCATE(Ez_phase_rec(xn+1,yn+1))
 !      ALLOCATE(Et_abs_rec(xn+1,yn+1))
 !      ALLOCATE(Emean_rec(xn+1,yn+1))
 ! ENDIF
END SUBROUTINE COM_RCVRGRID
!888888888888888888888888888888888
SUBROUTINE COM_FLOOR()
IMPLICIT NONE
CHARACTER(len=120) :: line1temp
REAL :: maxz,maxy,maxx,minz,miny,minx
INTEGER :: i,jj
maxz = point_mat(1,3)
minz = point_mat(1,3)
maxy = point_mat(1,2)
miny = point_mat(1,2)
maxx = point_mat(1,1)
minx = point_mat(1,1)
DO i=2,pnt_num
        IF (point_mat(i,3)>maxz) maxz=point_mat(i,3)
        IF (point_mat(i,3)<minz) minz=point_mat(i,3)
        IF (point_mat(i,2)>maxy) maxy=point_mat(i,2)
        IF (point_mat(i,2)<miny) miny=point_mat(i,2)
        IF (point_mat(i,1)>maxx) maxx=point_mat(i,1)
        IF (point_mat(i,1)<minx) minx=point_mat(i,1)
END DO

pnt_num = pnt_num + 1
line1temp = line1
WRITE(*,*) 'line1 = ', line1
WRITE(line1,*) 'floor1  ', minx , miny, minz
WRITE(*,*) 'line1 = ', line1
CALL COM_POINT(pnt_num)
WRITE(*,*)'the point is ', (point_mat(pnt_num,jj),jj=1,3)
pnt_num = pnt_num + 1
WRITE(line1,*) 'floor2  ', maxx, '  ', maxy, '  ',minz
WRITE(*,*) 'line1 = ', line1
CALL COM_POINT(pnt_num)
WRITE(*,*)'the point is ', (point_mat(pnt_num,jj),jj=1,3)
wlltp_num = wlltp_num + 1
WRITE(line1,*) 'WTFLOOR   ','SOLID   ', TRIM(line1temp)
WRITE(*,*) 'line1 = ', line1
CALL COM_WALL_TYPE(wlltp_num)
WRITE(line1,*) 'floor1   floor2   WTFLOOR'
wll_def_num = wll_def_num+1
CALL COM_DEFINE_WALL(wll_def_num,pnt_num)
line1 = line1temp
END SUBROUTINE COM_FLOOR
!8888888888888888888888888888888
SUBROUTINE COM_CEILING()
IMPLICIT NONE
CHARACTER(len=120) :: line1temp,line2temp
CHARACTER(len=10) :: tempchar
REAL :: maxz,maxy,maxx,minz,miny,minx
INTEGER :: i,nn
maxz = point_mat(1,3)
minz = point_mat(1,3)
maxy = point_mat(1,2)
miny = point_mat(1,2)
maxx = point_mat(1,1)
minx = point_mat(1,1)
DO i=2,pnt_num
        IF (point_mat(i,3)>maxz) maxz=point_mat(i,3)
        IF (point_mat(i,3)<minz) minz=point_mat(i,3)
        IF (point_mat(i,2)>maxy) maxy=point_mat(i,2)
        IF (point_mat(i,2)<miny) miny=point_mat(i,2)
        IF (point_mat(i,1)>maxx) maxx=point_mat(i,1)
        IF (point_mat(i,1)<minx) minx=point_mat(i,1)
END DO

line1temp = line1
line2temp = line1
CALL GET_WORD_var(tempchar,nn)
DO i=1,120
        line2temp(i:i) = ' '
END DO
line2temp(1:10) = tempchar
WRITE(*,*) 'line2(1:40) = ',line2temp(1:40) 
CALL GET_WORD_var(tempchar,nn)
line2temp(15:24) = tempchar
WRITE(*,*) 'line2(1:40) = ',line2temp(1:40) 
WRITE(*,*) 'TRIM(line1) = ', TRIM(line1)
READ(line1,*) maxz
WRITE(*,*) 'maxz = ',maxz
WRITE(line1,*) 'ptceiling1  ', minx , miny, maxz
WRITE(*,*) 'pnt_num = ', pnt_num
pnt_num = pnt_num + 1
CALL COM_POINT(pnt_num)
pnt_num = pnt_num + 1
WRITE(line1,*) 'ptceiling2  ', maxx, maxy, maxz
WRITE(*,*) 'pnt_num = ', pnt_num
CALL COM_POINT(pnt_num)
wlltp_num = wlltp_num + 1
WRITE(line1,*) 'WTCEILING   ','SOLID   ', TRIM(line2temp)
CALL COM_WALL_TYPE(wlltp_num)
WRITE(line1,*) 'ptceiling1   ptceiling2   WTCEILING'
wll_def_num = wll_def_num+1
CALL COM_DEFINE_WALL(wll_def_num,pnt_num)
line1 = line1temp
END SUBROUTINE COM_CEILING
!88888888888888888888888888888888
SUBROUTINE readinput(inputfilename)
!USE commandresults
!USE ALL_DATA
CHARACTER(len=12), INTENT(IN) :: inputfilename
REAL :: a
INTEGER:: int_var,i,mat_num,j,jj

 CHARACTER(len=10), DIMENSION(2,3) :: charmat
 CHARACTER(len=10) :: string1,string2
 CHARACTER(len=20) :: command_name
 command_name= '                    '
! 'MATERIAL            '
! 'WALL_TYPE           '
! 'POINT               '
! 'DEFINE_WALL         '
! 'FREQUENCY           '
! 'SOURCE_LOCATION     '
! 'SOURCE_POWER        '
! 'DIPOLE              '
! 'RCVRPOINT           '
! 'RCVRLINE            '
! 'RCVRGRID            '
! '                    '
! '                    '
! '                    '
! '                    '
mat_num=0
wlltp_num=0
pnt_num=0
wll_def_num=0
OPEN(UNIT=20, FILE=inputfilename, STATUS='OLD',ACTION='READ',IOSTAT=int_var)
READ(20,'(A)',IOSTAT=int_var) line1
!WRITE(*,*)'IOSATA= ',int_var
! reading the command name
wholefile: DO WHILE (int_var == 0)
       !Read first line
       command_name= '                    '
       ! ignoring comments. they start with "c"
       DO WHILE (line1(1:1)=='C' .OR. line1(1:1)=='c')
              READ(20,'(A)',IOSTAT=int_var) line1
                !WRITE(*,*)'inside while in main code IOSATA= ',int_var
                IF (int_var /= 0) EXIT
       END DO
       ! ignoring blanks
       delete_blanks: DO
              IF (line1(1:1) /= ' ') EXIT delete_blanks
              line1 = line1(2:)
       END DO delete_blanks
       ! start command_name until see a blank
       i = 1
       DO
              command_name(i:i)=line1(1:1)
              line1 = line1(2:)
              IF (line1(1:1) == ' ') EXIT
              i = i+1
       END DO
        !WRITE(*,*)'line one in main code is ',line1
        !WRITE(*,*)'my command name is :',command_name,'x'

       selecting_command: SELECT CASE ((command_name))
        CASE('MATERIAL            ')
                !WRITE(*,*)'Does it come  MATERIAL?'
              mat_num = mat_num+1
              CALL COM_MATERIAL(mat_num)
       CASE('WALL_TYPE           ')
              !WRITE(*,*)'Does it come  WALL_TYPE?'
              wlltp_num = wlltp_num+1
              CALL COM_WALL_TYPE(wlltp_num)
       CASE('POINT               ')
              !WRITE(*,*)'Does it come  POINT?'
              pnt_num = pnt_num +1
              CALL COM_POINT(pnt_num)
       CASE('DEFINE_WALL         ')
              !WRITE(*,*)'Does it come  DEFINE_WALL?'
              wll_def_num = wll_def_num+1
              CALL COM_DEFINE_WALL(wll_def_num,pnt_num)
       CASE('FREQUENCY           ')
              !WRITE(*,*)'Does it come  FREQUENCY?'
              CALL COM_FREQUENCY()
       CASE('SOURCE_LOCATION     ')
              CALL COM_SOURCE_LOCATION(pnt_num)
       CASE('SOURCE_POWER        ')
              CALL COM_SOURCE_POWER()
       CASE('DIPOLE              ')
              !WRITE(*,*)'Does it come  DIPOLE?'
              CALL COM_DIPOLE()
       CASE('RCVRPOINT           ')
              CALL COM_RCVRPOINT(line1)
       CASE('RCVRLINE            ')
              CALL COM_RCVRLINE(pnt_num)
       CASE('RCVRGRID            ')
              CALL COM_RCVRGRID(pnt_num)
       CASE('FLOOR               ')
               CALL COM_FLOOR()
       CASE('EILING              ')
               CALL COM_CEILING()
               
       CASE('                    ')
       END SELECT selecting_command
! if in subroutine calltest(vartest), vartest is CHARACTER(8)
! and if we call calltest('12345'); then stored in leftmost. good
       READ(20,'(A)',IOSTAT=int_var) line1
       !WRITE(*,*)'IOSATA= ',int_var
END DO wholefile
     n_mat = mat_num
     W = wll_def_num
     n_wt = wlltp_num
     ALLOCATE(Materials(n_mat,3))
     ALLOCATE(Walls(W,7))
     ALLOCATE(WTypes(n_wt,15))

     DO i=1,n_mat
       Materials(i,:)=material_mat(i,:)
     ENDDO
     
     DO i=1,n_wt
       WTypes(i,1)=wall_type_mat(i,1)
        DO j=2,INT(WTypes(i,1))+1
              WTypes(i,j)=wall_type_mat(i,1+(j-1-1)*2+2) !prone to error
       ENDDO
        DO j=(INT(WTypes(i,1))+2),(2*INT(WTypes(i,1))+1)
              !write(*,*)'j=',j,'(.)=',1+(j-1-1)*2+1,'wall_type_mat(.)=',wall_type_mat(i,1+(j-1-1)*2+1)
              jj=j-WTypes(i,1)
                !Write(*,*) 'jj =',jj
              WTypes(i,j)=wall_type_mat(i,1+(jj-1-1)*2+1) !prone to error
        ENDDO
     ENDDO
     
     DO i=1,wll_def_num
       Walls(i,1:3)=(/point_mat(define_wall_mat(i,1),:) /)! , point_mat(define_wall_mat(i,2),:), define_wall_types(i)/)
       Walls(i,4:6)=(/ point_mat(define_wall_mat(i,2),:) /)
       Walls(i,7) = define_wall_types(i)
     ENDDO
     DO i=1,n_mat
       !Write(*,*) WTypes(i,:)
     ENDDO
     !WRITE(*,*) 'frequency =',TX(5)
     !WRITE(*,*) 'TX location =',TX(1:3)
     !WRITE(*,*) 'power in mW =',TX(4)
     !WRITE(*,*) 'xrec =',xrec
     !WRITE(*,*) 'yrec =',yrec
     !WRITE(*,*) 'zrec =',zrec
END SUBROUTINE readinput
END MODULE commandresults

!888888888888888888888888888888888888888888
!888888888888888888888888888888888888888888
!888888888888888888888888888888888888888888
MODULE realloc_mod
CONTAINS
  FUNCTION reallocate_r2(p, n)               ! reallocate REAL
    REAL, POINTER, DIMENSION(:,:) :: p, reallocate_r2
    INTEGER, DIMENSION(2),intent(in) :: n
    INTEGER :: nold, ierr
    ALLOCATE(reallocate_r2(n(1),n(2)), STAT=ierr)
    IF(ierr /= 0) STOP "allocate error"
    IF(.NOT. ASSOCIATED(p)) RETURN
    
    
    DEALLOCATE(p) 
  END FUNCTION REALLOCATE_r2

  FUNCTION reallocate_r3(p, n)               ! reallocate REAL
    REAL, POINTER, DIMENSION(:,:,:) :: p, reallocate_r3
    INTEGER, DIMENSION(3),intent(in) :: n
    INTEGER :: nold, ierr
    ALLOCATE(reallocate_r3(n(1),n(2),n(3)), STAT=ierr)
    IF(ierr /= 0) STOP "allocate error"
    IF(.NOT. ASSOCIATED(p)) RETURN
    DEALLOCATE(p) 
  END FUNCTION REALLOCATE_r3

  FUNCTION reallocate_r1(p, n)               ! reallocate REAL
    REAL, POINTER, DIMENSION(:) :: p, reallocate_r1
    INTEGER, intent(in) :: n
    INTEGER :: nold, ierr
    ALLOCATE(reallocate_r1(n), STAT=ierr)
    IF(ierr /= 0) STOP "allocate error"
    IF(.NOT. ASSOCIATED(p)) RETURN
    DEALLOCATE(p) 
  END FUNCTION REALLOCATE_r1

  FUNCTION reallocate_c2(p, n)               ! reallocate REAL
    COMPLEX, POINTER, DIMENSION(:,:) :: p, reallocate_c2
    INTEGER, DIMENSION(2),intent(in) :: n
    INTEGER :: nold, ierr
    ALLOCATE(reallocate_c2(n(1),n(2)), STAT=ierr)
    IF(ierr /= 0) STOP "allocate error"
    IF(.NOT. ASSOCIATED(p)) RETURN
    DEALLOCATE(p) 
  END FUNCTION REALLOCATE_c2

END MODULE realloc_mod

!8888888888888888888888888888888888888888888888
!8888888888888888888888888888888888888888888888
!8888888888888888888888888888888888888888888888

