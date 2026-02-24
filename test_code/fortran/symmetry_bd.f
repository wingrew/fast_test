subroutine symmetry_bd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1),-ord+1:extc(2),-ord+1:extc(3)),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+1,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+1)*SoA(3)
   enddo

end subroutine symmetry_bd