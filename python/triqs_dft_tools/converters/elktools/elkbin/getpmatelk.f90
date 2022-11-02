
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getpmatelk(ik,nstsv,vkl,pmat)
!use modmain
implicit none
! arguments
!real(8), intent(in) :: vpl(3)
integer, intent(in) :: ik !Need to read this in for the interface
integer, intent(in) :: nstsv !Need to read this in for the interface
real(8), intent(in) :: vkl(3) !TRIQS uses reduced kpts set
complex(8), intent(out) :: pmat(nstsv,nstsv,3)
! local variables
!integer isym,ik,ist,jst
integer recl,nstsv_,i
real(8) vkl_(3),t1
!real(8) vkl_(3),sc(3,3),t1
!complex(8) v1(3),v2(3)

!adnj - set up tolerance for lattice vectors, although this is an input in Elk,
! it is not advised to change this in Elk. Therefore, it should be fine to set 
!it as a constant here.
real(8) epslat
epslat=1.d-6

! find the equivalent k-point number and symmetry which rotates vkl to vpl
! adnj edit - not needed for interface: TRIQS uses reduced kpts set.
!call findkpt(vpl,isym,ik)
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
!t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
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

! adnj - not needed as TRIQS uses reduced kpts set
! if p = k then return
!t1=abs(vpl(1)-vkl(1,ik))+abs(vpl(2)-vkl(2,ik))+abs(vpl(3)-vkl(3,ik))
!if (t1.lt.epslat) return
! rotate the matrix elements from the reduced to non-reduced k-point
!sc(:,:)=symlatc(:,:,lsplsymc(isym))
!do ist=1,nstsv
!  do jst=1,nstsv
!    v1(:)=pmat(ist,jst,:)
!    call rz3mv(sc,v1,v2)
!    pmat(ist,jst,:)=v2(:)
!  end do
!end do
return

!This subroutine is not needed for TRIQS interface
!contains

!subroutine rz3mv(a,x,y)
!implicit none
!real(8), intent(in) :: a(3,3)
!complex(8), intent(in) :: x(3)
!complex(8), intent(out) :: y(3)
!y(1)=a(1,1)*x(1)+a(1,2)*x(2)+a(1,3)*x(3)
!y(2)=a(2,1)*x(1)+a(2,2)*x(2)+a(2,3)*x(3)
!y(3)=a(3,1)*x(1)+a(3,2)*x(2)+a(3,3)*x(3)
!return
!end subroutine

end subroutine

