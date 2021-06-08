module smd_io
   use mctc_env, only: wp
   use mctc_io_symbols, only: to_number
   private

   public :: read_smd

   integer, parameter :: max_elem=94

   interface read_smd
      module procedure read_solvent_properties
      module procedure read_parameters_h2o
      module procedure read_parameters_ot
   end interface read_smd

contains

   subroutine read_parameters_ot(filename,zk,zkk,rzkk,drzkk,nc3,rnc3,drnc3,&
         &sg,sr,sp,sb)
      !> Parameter filename for solvent Parameters
      character(len=*), intent(in) :: filename
      !> Read Solvent Parameters
      real(wp), intent(out) :: zk(3,max_elem),zkk(3,max_elem,max_elem),&
         &rzkk(max_elem,max_elem),drzkk(max_elem,max_elem),&
         &nc3,rnc3,drnc3,sg,sr,sp,sb

      !> read line
      character(len=100) :: line, lines(max_elem) 
      !> read Symbol
      character(len=2) :: symbol, dumsymb
      integer, allocatable :: symbol_id(:)
      real(wp), allocatable :: zkk_temp(:,:)
      real(wp) :: par
      !> Sanity Check
      logical :: ex
      !> Laufvariable
      integer :: i, j
      !> Zählvariable 
      integer :: elecount
      !> IO_error
      integer :: io_error

      INQUIRE(file=filename,exist=ex)
      if (.not. ex) then
         write(*,*) "Something went wrong, H2O Parameter file was specified, but does not exist."
         error stop
      end if
      
      !> Initialize not specified Params as zero
      zk=0.0_wp
      zkk=0.0_wp
      rzkk=0.0_wp
      drzkk=0.0_wp

      !> Open Parameter File
      open(unit=1,file=filename)

      read(1,*) ! Skip first line

      !> Read the specified Zk Values
      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#Zk2") exit
            read(line,*) symbol, par
            zk(1,to_number(symbol))=par
         else
            write(*,*) "Error while reading Zk Values in the H2O Parameter file."
            error stop
         end if
      end do

      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#Zk3") exit
            read(line,*) symbol, par
            zk(2,to_number(symbol))=par
         else
            write(*,*) "Error while reading Zk2 Values in the H2O Parameter file."
            error stop
         end if
      end do

      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#Zkk Matrix") exit
            read(line,*) symbol, par
            zk(3,to_number(symbol))=par
         else
            write(*,*) "Error while reading Zk3 Values in the H2O Parameter file."
            error stop
         end if
      end do

      !> Count Elements of the Zkk1 Matrix      

      elecount=0
      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#Zkk2 Matrix") exit
            elecount=elecount+1
            lines(elecount)=line
         else
            write(*,*) "Error while determining Elements of the Zkk Matrix."
            error stop
         end if
      end do

      !> Appropriately allocating Symbol ID Matrix
      
      allocate(symbol_id(elecount))
      allocate(zkk_temp(elecount,elecount))

      do i=1,elecount
         read(lines(i),*) symbol, zkk_temp(i,:)
         symbol_id(i)=to_number(symbol)
      end do

      do i=1,elecount
         do j=1,elecount
            zkk(1,symbol_id(i),symbol_id(j))=zkk_temp(i,j)
         end do
      end do

      deallocate(symbol_id, zkk_temp)

      !> Count Elements of the Zkk Matrix

      elecount=0
      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#Zkk3 Matrix") exit
            elecount=elecount+1
            lines(elecount)=line
         else
            write(*,*) "Error while determining Elements of the Zkk2 Matrix."
            error stop
         end if
      end do

      !> Appropriately allocating Symbol ID Matrix
      
      allocate(symbol_id(elecount))
      allocate(zkk_temp(elecount,elecount))

      do i=1,elecount
         read(lines(i),*) symbol, zkk_temp(i,:)
         symbol_id(i)=to_number(symbol)
      end do

      do i=1,elecount
         do j=1,elecount
            zkk(2,symbol_id(i),symbol_id(j))=zkk_temp(i,j)
         end do
      end do

      deallocate(symbol_id, zkk_temp)

      !> Count Elements of the Zkk3 Matrix

      elecount=0
      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#rzkk Matrix") exit
            elecount=elecount+1
            lines(elecount)=line
         else
            write(*,*) "Error while determining Elements of the Zkk3 Matrix."
            error stop
         end if
      end do

      !> Appropriately allocating Symbol ID Matrix
      
      allocate(symbol_id(elecount))
      allocate(zkk_temp(elecount,elecount))

      do i=1,elecount
         read(lines(i),*) symbol, zkk_temp(i,:)
         symbol_id(i)=to_number(symbol)
      end do

      do i=1,elecount
         do j=1,elecount
            zkk(3,symbol_id(i),symbol_id(j))=zkk_temp(i,j)
         end do
      end do

      deallocate(symbol_id, zkk_temp)

      !> Count Elements of the rzkk Matrix

      elecount=0
      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#drzkk Matrix") exit
            elecount=elecount+1
            lines(elecount)=line
         else
            write(*,*) "Error while determining Elements of the rzkk Matrix."
            error stop
         end if
      end do

      !> Appropriately allocating Symbol ID Matrix
      
      allocate(symbol_id(elecount))
      allocate(zkk_temp(elecount,elecount))

      do i=1,elecount
         read(lines(i),*) symbol, zkk_temp(i,:)
         symbol_id(i)=to_number(symbol)
      end do

      do i=1,elecount
         do j=1,elecount
            rzkk(symbol_id(i),symbol_id(j))=zkk_temp(i,j)
         end do
      end do

      deallocate(symbol_id, zkk_temp)

      !> Count Elements of the drzkk Matrix

      elecount=0
      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#NC3") exit
            elecount=elecount+1
            lines(elecount)=line
         else
            write(*,*) "Error while determining Elements of the drzkk Matrix."
            error stop
         end if
      end do

      !> Appropriately allocating Symbol ID Matrix
      
      allocate(symbol_id(elecount))
      allocate(zkk_temp(elecount,elecount))

      do i=1,elecount
         read(lines(i),*) symbol, zkk_temp(i,:)
         symbol_id(i)=to_number(symbol)
      end do

      do i=1,elecount
         do j=1,elecount
            drzkk(symbol_id(i),symbol_id(j))=zkk_temp(i,j)
         end do
      end do

      deallocate(symbol_id, zkk_temp)

      read(1,*) nc3
      read(1,*) rnc3
      read(1,*) drnc3
      read(1,*) !Solvent Line
      read(1,*) sg
      read(1,*) sr
      read(1,*) sp
      read(1,*) sb
   end subroutine read_parameters_ot


   subroutine read_parameters_h2o(filename,zk,zkk,rzkk,drzkk,nc3,rnc3,drnc3)
      !> Parameter filename for solvent Parameters
      character(len=*), intent(in) :: filename
      !> Read Solvent Parameters
      real(wp), intent(out) :: zk(max_elem),zkk(max_elem,max_elem),rzkk(max_elem,max_elem)&
         &,drzkk(max_elem,max_elem),nc3,rnc3,drnc3

      !> read line
      character(len=100) :: line, lines(max_elem) 
      !> read Symbol
      character(len=2) :: symbol, dumsymb
      integer, allocatable :: symbol_id(:)
      real(wp), allocatable :: zkk_temp(:,:)
      real(wp) :: par
      !> Sanity Check
      logical :: ex
      !> Laufvariable
      integer :: i, j
      !> Zählvariable 
      integer :: elecount
      !> IO_error
      integer :: io_error

      INQUIRE(file=filename,exist=ex)
      if (.not. ex) then
         write(*,*) "Something went wrong, H2O Parameter file was specified, but does not exist."
         error stop
      end if
      
      !> Initialize not specified Params as zero
      zk=0.0_wp
      zkk=0.0_wp
      rzkk=0.0_wp
      drzkk=0.0_wp

      !> Open Parameter File
      open(unit=1,file=filename)

      read(1,*) ! Skip first line

      !> Read the specified Zk Values
      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#Zkk Matrix") exit
            read(line,*) symbol, par
            zk(to_number(symbol))=par
         else
            write(*,*) "Error while reading Zk Values in the H2O Parameter file."
            error stop
         end if
      end do

      !> Count Elements of the Zkk Matrix      

      elecount=0
      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#rzkk Matrix") exit
            elecount=elecount+1
            lines(elecount)=line
         else
            write(*,*) "Error while determining Elements of the Zkk Matrix."
            error stop
         end if
      end do

      !> Appropriately allocating Symbol ID Matrix
      
      allocate(symbol_id(elecount))
      allocate(zkk_temp(elecount,elecount))

      do i=1,elecount
         read(lines(i),*) symbol, zkk_temp(i,:)
         symbol_id(i)=to_number(symbol)
      end do

      do i=1,elecount
         do j=1,elecount
            zkk(symbol_id(i),symbol_id(j))=zkk_temp(i,j)
         end do
      end do

      deallocate(symbol_id, zkk_temp)

      !> Count Elements of the rzkk Matrix

      elecount=0
      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#drzkk Matrix") exit
            elecount=elecount+1
            lines(elecount)=line
         else
            write(*,*) "Error while determining Elements of the rzkk Matrix."
            error stop
         end if
      end do

      !> Appropriately allocating Symbol ID Matrix
      
      allocate(symbol_id(elecount))
      allocate(zkk_temp(elecount,elecount))

      do i=1,elecount
         read(lines(i),*) symbol, zkk_temp(i,:)
         symbol_id(i)=to_number(symbol)
      end do

      do i=1,elecount
         do j=1,elecount
            rzkk(symbol_id(i),symbol_id(j))=zkk_temp(i,j)
         end do
      end do

      deallocate(symbol_id, zkk_temp)

      !> Count Elements of the drzkk Matrix

      elecount=0
      do while (.true.)
         read(1,'(A)',iostat=io_error) line
         if (io_error .EQ. 0) then
            if (trim(line) .EQ. "#NC3") exit
            elecount=elecount+1
            lines(elecount)=line
         else
            write(*,*) "Error while determining Elements of the drzkk Matrix."
            error stop
         end if
      end do

      !> Appropriately allocating Symbol ID Matrix
      
      allocate(symbol_id(elecount))
      allocate(zkk_temp(elecount,elecount))

      do i=1,elecount
         read(lines(i),*) symbol, zkk_temp(i,:)
         symbol_id(i)=to_number(symbol)
      end do

      do i=1,elecount
         do j=1,elecount
            drzkk(symbol_id(i),symbol_id(j))=zkk_temp(i,j)
         end do
      end do

      deallocate(symbol_id, zkk_temp)

      read(1,*) nc3
      read(1,*) rnc3
      read(1,*) drnc3
   end subroutine read_parameters_h2o

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

      close(1)

   end subroutine read_solvent_properties


end module smd_io
