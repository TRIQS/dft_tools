! Regression test for the dmftproj spin-orbit symmetry construction.
!
! The spinor symmetry matrices written by dmftproj (built from dmat/d_matrix and
! the orbital time-reversal operator timeinv_op) must form a valid (anti)unitary
! representation of the magnetic point group. The dft_tools symmetrizer then
! averages M_sym = (1/N) sum_g  D_g (M or conj M) D_g^dag and must be an
! idempotent projector, so symmetrizing an already-symmetric matrix leaves its
! eigenvalues unchanged. This drives the real dmat and timeinv_op on a D4
! magnetic group (z-rotations plus time-reversed in-plane C2 axes) and checks
! idempotency for l=1 and l=2. See TRIQS/dft_tools#148.
program test_spinrot_symmetry
  implicit none
  integer, parameter :: dp = 8
  integer, parameter :: ng = 8
  real(dp), parameter :: tol = 1.0d-8
  ! D4 magnetic group, z-y-z Euler angles a,b,g, 3D determinant, time-reversal flag, phase
  real(dp) :: ad(ng), bd(ng), gd(ng), dt(ng), ph(ng)
  integer  :: ti(ng)
  logical  :: ok
  data ad /0d0, 1.5707963267948966d0, 3.1415926535897931d0, -1.5707963267948966d0, &
           3.1415926535897931d0, 0d0, -1.5707963267948966d0, 1.5707963267948966d0/
  data bd /0d0,0d0,0d0,0d0, 3.1415926535897931d0,3.1415926535897931d0, &
           3.1415926535897931d0,3.1415926535897931d0/
  data gd /0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/
  data dt /1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/
  data ti /0,0,0,0,1,1,1,1/
  data ph /0d0,1.5707963267948966d0,3.1415926535897931d0,-1.5707963267948966d0, &
           -3.1415926535897931d0,0d0,1.5707963267948966d0,-1.5707963267948966d0/

  ok = .true.
  call check(1, ok)
  call check(2, ok)
  if (.not. ok) then
    write(*,'(A)') 'FAILED: symmetrizer is not idempotent for the production construction'
    error stop 1
  end if
  write(*,'(A)') 'PASSED'

contains

  subroutine check(l, ok)
    integer, intent(in) :: l
    logical, intent(inout) :: ok
    real(dp) :: idem
    idem = idempotency(l, 1)        ! production sign: ephase = exp(+i*phase/2)
    write(*,'(A,I2,A,ES10.3)') 'l=', l, '  idempotency |S^2-S| = ', idem
    if (idem > tol) ok = .false.
  end subroutine

  function idempotency(l, sgn) result(res)
    integer, intent(in) :: l, sgn
    real(dp) :: res
    integer :: d
    complex(dp), allocatable :: M0(:,:), Ms(:,:), Mss(:,:)
    d = 2*(2*l+1)
    allocate(M0(d,d), Ms(d,d), Mss(d,d))
    call herm_seed(d, M0)
    call symmetrize(l, sgn, M0, Ms)
    call symmetrize(l, sgn, Ms, Mss)
    res = maxval(abs(Mss - Ms))
    deallocate(M0, Ms, Mss)
  end function

  subroutine symmetrize(l, sgn, M, out)
    integer, intent(in) :: l, sgn
    complex(dp), intent(in)  :: M(:,:)
    complex(dp), intent(out) :: out(:,:)
    integer :: ig, d
    complex(dp), allocatable :: mat(:,:)
    d = 2*(2*l+1)
    allocate(mat(d,d))
    out = (0d0,0d0)
    do ig = 1, ng
      call build_mat(l, ig, sgn, mat)
      if (ti(ig) == 1) then
        out = out + matmul(matmul(mat, conjg(M)), conjg(transpose(mat))) / real(ng, dp)
      else
        out = out + matmul(matmul(mat, M), conjg(transpose(mat))) / real(ng, dp)
      end if
    end do
    deallocate(mat)
  end subroutine

  subroutine build_mat(l, ig, sgn, mat)
    integer, intent(in) :: l, ig, sgn
    complex(dp), intent(out) :: mat(:,:)
    integer :: n
    complex(dp), allocatable :: rotl(:,:)
    complex(dp) :: eph
    n = 2*l + 1
    allocate(rotl(n,n))
    call dmat(l, ad(ig), bd(ig), gd(ig), dt(ig), rotl, n)
    if (ti(ig) == 1) call timeinv_op(rotl, n, l, 0)
    eph = exp(cmplx(0d0, real(sgn,dp)*ph(ig)/2d0, dp))
    mat = (0d0,0d0)
    mat(1:n, 1:n)         = eph * rotl
    mat(n+1:2*n, n+1:2*n) = conjg(eph) * rotl
    deallocate(rotl)
  end subroutine

  subroutine herm_seed(d, M)
    integer, intent(in) :: d
    complex(dp), intent(out) :: M(d,d)
    integer :: i, j
    do i = 1, d
      do j = 1, d
        M(i,j) = cmplx(sin(0.7d0*i + 1.3d0*j), cos(0.4d0*i - 0.9d0*j), dp)
      end do
    end do
    M = 0.5d0*(M + conjg(transpose(M)))
  end subroutine

end program
