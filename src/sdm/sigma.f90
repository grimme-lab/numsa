module sdm_sigma
   use sdm_init, only: sdm_param
   use mctc_env, only: wp
   use mctc_io_symbols, only: to_number
   use mctc_io_convert, only: autoaa
   implicit none
   private
   public :: calc_surft, sdm_surft

   integer, parameter :: max_elem=94

   interface calc_surft
      module procedure :: calc_surft_h2o
      !module procedure :: calc_surft
   end interface

   type :: sdm_info
      !> Number of Atoms in the Solute
      integer :: nat
      !> Identifier of Atom i
      integer, allocatable :: Z(:)
   end type sdm_info

   type :: sdm_surft
      !> Atom depended surface tension
      real(wp) :: sk(max_elem)
      !> Solvent depended molecular surface tension
      real(wp) :: sm
   end type

contains

   subroutine calc_surft_h2o(xyz,species,ident,param,surft)
      !> Atomic Coordinates of the Solute in bohr [3,nat]
      real(wp), intent(in) :: xyz(:,:)
      !> Unique Chemical Species in the Solute [nat]
      integer, intent(in) :: species(:)
      !> Identifier (Element Symbol) of the unique chemical species in the solute [nat]
      character(len=*), intent(in) :: ident(:)
      !> SDM Parameters
      type(sdm_param), intent(in) :: param
      !> SDM Surft Output
      type(sdm_surft), intent(out) :: surft

      !> Local Variables
      !> Needed Info for sdm
      type(sdm_info) :: self
      !> Laufvariabeln
      integer :: Z, i, j, k
      !> Temporary sigma saving
      real(wp) :: s_temp1, s_temp2, s_temp3, s_temp4, nc_temp
      
      Call init_info(xyz,species,ident,self)

      s_temp1=0.0_wp
      s_temp2=0.0_wp
      s_temp3=0.0_wp
      s_temp4=0.0_wp
      nc_temp=0.0_wp

      surft%sm=0.0_wp !For H2O, sm is zero
 
      ! Setting up s_k Parameters now
      do Z=1,max_elem
         select case(Z)
            case (1) !H
               do i=1,self%nat
                  if (self%Z(i) .EQ. 1) then
                     do j=1, self%nat
                        select case(self%Z(j))
                           case (6) !H,C
                              s_temp1=s_temp1+T(xyz(:,i),xyz(:,j),param%rzkk(1,6),param%drzkk(1,6))
                           case (8) !H,O
                              s_temp2=s_temp2+T(xyz(:,i),xyz(:,j),param%rzkk(1,8),param%drzkk(1,8))
                           case default
                              cycle
                        end select
                     end do
                  end if
               end do
               surft%sk(1)=param%zk(1,1)+param%zkk(1,1,6)*s_temp1+param%zkk(1,1,8)*s_temp2
               s_temp1=0.0_wp
               s_temp2=0.0_wp
            case (6) !C
               do i=1,self%nat
                  if (self%Z(i) .EQ. 6) then
                     do j=1, self%nat
                        select case(self%Z(j))
                           case (6) !C,C
                              if (i .NE. j) s_temp1=s_temp1+T(xyz(:,i),xyz(:,j),param%rzkk(6,6),param%drzkk(6,6))
                           case (7) !C,N
                              s_temp2=s_temp2+T(xyz(:,i),xyz(:,j),param%rzkk(6,7),param%drzkk(6,7))
                           case default
                              cycle
                        end select
                     end do
                  end if
               end do
               surft%sk(6)=param%zk(1,6)+param%zkk(1,6,6)*s_temp1+param%zkk(1,6,7)*(s_temp2**2)
               s_temp1=0.0_wp
               s_temp2=0.0_wp
            case(7) !N
               do i=1,self%nat
                  if (self%Z(i) .EQ. 7) then
                     do j=1, self%nat
                        select case(self%Z(j))
                           case (6) !N,C
                              s_temp1=s_temp1+T(xyz(:,i),xyz(:,j),param%rzkk(6,6),param%drzkk(6,6))
                              do k=1, self%nat ! For N, all interactions between C and all other atoms play a role
                                 if ((k .NE. j) .AND. (k .NE. i)) then
                                    nc_temp=nc_temp+T(xyz(:,j),xyz(:,k),param%rzkk(6,self%Z(k)),param%drzkk(6,self%Z(k)))
                                 end if
                              end do
                              s_temp1=s_temp1*(nc_temp**2)
                              nc_temp=0.0_wp
                              s_temp2=s_temp2+T(xyz(:,i),xyz(:,j),param%rnc3,param%drnc3) !Additional nc3 parametric correction
                           case default
                              cycle
                        end select
                     end do
                  end if
               end do
               surft%sk(7)=param%zk(1,7)+param%zkk(1,7,6)*(s_temp1**1.3_wp)+param%nc3*s_temp2
               s_temp1=0.0_wp
               s_temp2=0.0_wp
            case(8) !O
               do i=1,self%nat
                  if (self%Z(i) .EQ. 8) then
                     do j=1, self%nat
                        select case(self%Z(j))
                           case (6) !O,C
                              s_temp1=s_temp1+T(xyz(:,i),xyz(:,j),param%rzkk(8,6),param%drzkk(8,6))
                           case (7) !O,N
                              s_temp2=s_temp2+T(xyz(:,i),xyz(:,j),param%rzkk(8,7),param%drzkk(8,7))
                           case (8) !O,O
                              if (i .NE. j) s_temp3=s_temp3+T(xyz(:,i),xyz(:,j),param%rzkk(8,8),param%drzkk(8,8))
                           case (15) !O,P
                              s_temp4=s_temp4+T(xyz(:,i),xyz(:,j),param%rzkk(8,15),param%drzkk(8,15))
                           case default
                              cycle
                        end select
                     end do
                  end if
               end do
               surft%sk(8)=param%zk(1,8)+param%zkk(1,8,6)*s_temp1+param%zkk(1,8,7)&
                  &*s_temp2+param%zkk(1,8,8)*s_temp3+param%zkk(1,8,15)*s_temp4
               s_temp1=0.0_wp
               s_temp2=0.0_wp
               s_temp3=0.0_wp
               s_temp4=0.0_wp
            case(9,14,16,17,35) !F, Si, S, Cl, Br
               surft%sk(Z)=param%zk(1,Z)
            case default
               surft%sk(Z)=0.0_wp !The default case for the surface tension is zero
         end select
      end do

   end subroutine calc_surft_h2o

   !subroutine calc_surft(xyz,n,alpha,beta,surft,arom,fclbr,sdm,symbol,zk_in,zkk_in,nc3_in)
      !!>Atom coordinates of the Solute
      !real(wp), intent(in) :: xyz(:,:)
      !!> Refraction index, Abrahams hydrogen bond accidity and Abrahams hydrogen bond basicity of the Solvent
      !real(wp), intent(in) :: n, alpha, beta
      !!> macroscopic surface tension, fraction of aromatic atoms and fraction of f,cl and br of the solvent
      !real(wp), intent(in) :: surft,arom,fclbr
      !!> Self defined Parameter for the SDM
      !real(wp), intent(in), optional :: zk_in(:), zkk_in(:,:),nc3_in
      !!> Symbols for self defined Parameters
      !character(len=*), intent(in), optional :: symbol
      !!> SDM Parameters
      !type(sdm_surft), intent(out) :: sdm

      !Call init_default(.false.)


   !end subroutine init_sdm

   function T(xyz1,xyz2,rzkk,drzkk)
      !>Atom coordinates xyz [3]
      real(wp), intent(in) :: xyz1(3), xyz2(3)
      !>rzkk and drzkk Parameters for the two Elements
      real(wp), intent(in) :: rzkk, drzkk
      !> Output of the T switching Function
      real(wp) :: T

      !>Distance between two Atoms
      real(wp) :: R

      !>Temporary for denominator
      real(wp) :: denom

      R=sqrt((xyz1(1)-xyz2(1))**2+(xyz1(2)-xyz2(2))**2+(xyz1(3)-xyz2(2))**2)**autoaa

      if (R .LT. (rzkk+drzkk)) then
         denom=R-drzkk-rzkk
         T=exp(drzkk/denom)
      else
         T=0.0_wp
      end if

   end function T

   subroutine init_info(xyz,species,ident,self)
      !> Atomic Coordinates of the Solute [3,nat]
      real(wp), intent(in) :: xyz(:,:)
      !> Unique Chemical Species in the Solute [nat]
      integer, intent(in) :: species(:)
      !> Identifier (Element Symbol) of the unique chemical species in the solute [nat]
      character(len=*), intent(in) :: ident(:)
      !> Output Info
      type(sdm_info) :: self

      !> Laufvariable
      integer :: elem
      self%nat=size(species)
      allocate(self%Z(self%nat))

      do elem=1,self%nat
         self%Z(elem)=to_number(ident(species(elem)))
      end do
   end subroutine init_info 
      


end module sdm_sigma
