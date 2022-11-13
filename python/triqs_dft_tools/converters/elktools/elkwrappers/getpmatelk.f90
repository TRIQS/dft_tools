
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.

! version 6.2.8 file modified by A. D. N. James for interface with TRIQS

subroutine getpmatelk(ik,nstsv,vkl,pmat)
!use modmain
implicit none
! arguments
integer, intent(in) :: ik !Need to read this in for the interface
integer, intent(in) :: nstsv !Need to read this in for the interface
real(8), intent(in) :: vkl(3) !TRIQS uses reduced kpts set
complex(8), intent(out) :: pmat(nstsv,nstsv,3)
! local variables
integer recl,nstsv_,i
real(8) vkl_(3),t1

!adnj - set up tolerance for lattice vectors, although this is an input in Elk,
! it is not advised to change this in Elk. Therefore, it should be fine to set 
!it as a constant here.
real(8) epslat
epslat=1.d-6

! find the record length
inquire(iolength=recl) vkl_,nstsv_,pmat
!$OMP CRITICAL(u150)
do i=1,2
  open(150,file='PMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl,err=10)
  read(150,rec=ik,err=10) vkl_,nstsv_,pmat
  exit
10 continue
  if (i.eq.2) then
    write(*,*)
    write(*,'("Error(getpmat): unable to read from PMAT.OUT")')
    write(*,*)
    stop
  end if
  close(150)
end do
!$OMP END CRITICAL(u150)
!adnj edit - updated for vkl array from TRIQS
t1=abs(vkl(1)-vkl_(1))+abs(vkl(2)-vkl_(2))+abs(vkl(3)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getpmat): differing vectors for k-point ",I8)') ik
  !write(*,'(" current  : ",3G18.10)') vkl(:,ik)
  write(*,'(" current  : ",3G18.10)') vkl(:)
  write(*,'(" PMAT.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getpmat): differing nstsv for k-point ",I8)') ik
  write(*,'(" current  : ",I8)') nstsv
  write(*,'(" PMAT.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if

return

end subroutine

