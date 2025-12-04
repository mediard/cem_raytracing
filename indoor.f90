PROGRAM indoor
USE raytracingmod
!USE sabinemod
IMPLICIT NONE
CALL raytracinggMod('myroommm.go3')
! It is in raytracing that everything if initialized.
! So Sabine SUBR. should be always called after raytracing.
!WRITE(*,*) 'xn = ', xn, '  yn = ', yn
!READ(*,*)
!CALL Sabine('flr_plan.go3')
END PROGRAM indoor
