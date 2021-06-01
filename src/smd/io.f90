module smd_io
   use mctc_env, only: wp
   private

   public :: read_smd

   interface read_smd
      module procedure read_solvent_properties
   end interface read_smd

contains

   subroutine read_solvent_properties(filename,n,alpha,beta,msurft,arom,fclbr)
      !> Propertyfile for the Solvent Properties
      character(len=*), intent(in) :: filename
      !> Output Properties
      real(wp), intent(out) :: n, alpha, beta, msurft, arom, fclbr

      !> Sanity Check
      logical :: ex
      !> Laufvariable
      integer :: i
      !> IO_Error
      integer :: io_error

      INQUIRE(file=filename,exist=ex)
      if (.not. ex) then
         write(*,*) "Something went wrong, property file does not seem to exist."
         error stop
      end if
      open(unit=1,file=filename)

      read(1,*,iostat=io_error) n, alpha, beta, msurft, arom, fclbr

      if (io_error .ne. 0) then
         write(*,*) "Something went wrong, property file may not have the right format."
         error stop
      end if

   end subroutine read_solvent_properties


end module smd_io
