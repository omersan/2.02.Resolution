!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Wavenumber plots
!     
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Aug. 27, 2015
!-----------------------------------------------------------------------------!

program wavenumbers
implicit none
integer::j,np
real*8 ::pi,hk,c2,c4,c6,p4,p6

pi = 4.0d0*datan(1.0d0)

! Writing modified wavenumbers for difference schemes
np = 2000 !use 2000 points to plot curve
open(12, file="wanenumbers.plt")
write(12,*)'variables ="hk","c2","c4","c6","Pade4","Pade6","exact"'
	do j=0,np
		hk = dfloat(j)*pi/dfloat(np)
		c2 = dsin(hk)
        c4 = (8.0d0*dsin(hk) - dsin(2.0d0*hk))/6.0d0
        c6 = (45.0d0*dsin(hk) - 9.0d0*dsin(2.0d0*hk) + dsin(3.0d0*hk) )/30.0d0
        p4 = 3.0d0*dsin(hk)/(2.0d0 + dcos(hk))
        p6 = (14.0d0/3.0d0*dsin(hk) + 1.0d0/6.0d0*dsin(2.0d0*hk))/(3.0d0 + 2.0d0*dcos(hk))
        
		write(12,*) hk,c2,c4,c6,p4,p6,hk
        
	end do
close(12)
     
end


