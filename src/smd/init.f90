module smd_init

!--------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------!
!-----               This Module initiates the smd Parameters               -----!
!-----              Own Parameters need to be added as follows:             -----!
!-----                                                                      -----!
!-----                             Symbol zk                                -----!
!-----                             Symbol zk                                -----!
!-----                             Symbol zk                                -----!
!-----                             #Zkk Matrix                              -----!
!----- Do not add this line!            Symbol Symbol Symbol                -----!
!-----                            Symbol zkk    zkk    zkk                  -----!
!-----                            Symbol zkk    zkk    zkk                  -----!
!-----                            Symbol zkk    zkk    zkk                  -----! 
!-----                            #NC3=nc3                                  -----!
!--------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------!
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   use mctc_io_convert, only : aatoau, autoaa
   use smd_io, only : read_smd
   implicit none
   private
   public :: init_smd, smd_param

   integer, parameter :: max_elem = 94
   
   type :: smd_param
      !> Element dependent Parameters
      real(wp) :: zk(max_elem)
      !> Pair/Solvent dependent Parameters
      real(wp) :: zkk(max_elem,max_elem)
      !> Pair dependent Parameters
      real(wp) :: rzkk(max_elem,max_elem)
      real(wp) :: drzkk(max_elem,max_elem)
      !> Additional nc3 Parameters
      real(wp) :: nc3
      real(wp) :: rnc3
      real(wp) :: drnc3
      !> Additional Parameters for Solvent Surface Tension
      real(wp) :: s_m !Macroscopic Solvent Surface Tension
   end type smd_param

 
   real(wp) :: ref_zk_h2o(max_elem)
   
   real(wp) :: ref_zkk_h2o(max_elem,max_elem)

   real(wp) :: ref_zk(3,max_elem)

   real(wp) :: ref_zkk(3,max_elem,max_elem)

   real(wp) :: ref_rzkk(max_elem,max_elem), ref_drzkk(max_elem, max_elem)

   real(wp) :: ref_nc3, ref_rnc3, ref_drnc3

   real(wp), parameter :: ref_sg=0.35_wp

   real(wp), parameter :: ref_sr2=-4.19_wp

   real(wp), parameter :: ref_sp2=-6.68_wp

   real(wp), parameter :: ref_sb2=0.00_wp

   !interface init_smd
      !module procedure :: init_smd_h2o
      !module procedure :: init_smd
   !end interface

contains

   subroutine init_smd(param,solvent)
      !>Solvent Name
      character(len=*), intent(in) :: solvent
      !>Output Parameters
      type(smd_param), intent(out) :: param

      !>Local Variables
      !>Are there Solvent Properties for the Solvent?
      logical :: ex
      !> Variables for Solvent Properties
      real(wp) :: n,alpha,beta,msurft,arom,fclbr

      ex=.false.

      select case(solvent)
         case('h2o', 'water')
            Call init_smd_h2o(param)
         case('methanol','ch4')
            Call init_smd_ot(param,1.3314_wp,0.43_wp,0.47_wp,32.3847_wp,0.0_wp,0.0_wp)
         case('dmso','DMSO')
            Call init_smd_ot(param,1.4772_wp,0.0_wp,0.88_wp,62.6680_wp,0.0_wp,0.0_wp)
         case('acetonitrile')
            Call init_smd_ot(param,1.3421_wp,0.07_wp,0.32_wp,28.4000_wp,0.0_wp,0.0_wp)
         case default
            INQUIRE(file=solvent//".prop",exist=ex)
            if (ex) then
               write(*,*) 'Reading self defined solvent properties for ', solvent,'.'
               Call read_smd(solvent//".prop",n,alpha,beta,msurft,arom,fclbr)
               Call init_smd_ot(param,n,alpha,beta,msurft,arom,fclbr)
            else
               write(*,*) 'No Solvent Properties for ', solvent, '.'
               stop
            end if
      end select

   end subroutine init_smd

   subroutine init_smd_h2o(param)
      !> smd Parameters
      type(smd_param), intent(out) :: param

      !> Local Variables
      !> Working Parameters
      real(wp) :: zk(max_elem), nc3, rnc3, drnc3
      real(wp), dimension(max_elem,max_elem) :: zkk, rzkk, drzkk

      !> Is there a self defined Parameter file?
      logical :: ex

      Inquire(file='smd_h2o.param',exist=ex)

      if (ex) then
         write(*,*) 'Self defined smd Parameters not implemented yet.'
         stop
      else
         write(*,*) "Using default smd Parameters."
         Call init_default(.TRUE.)
         param%zk(:)=ref_zk_h2o
         param%zkk(:,:)=ref_zkk_h2o
         param%rzkk=ref_rzkk
         param%drzkk=ref_drzkk
         param%nc3=ref_nc3
         param%rnc3=ref_rnc3
         param%drnc3=ref_drnc3
         param%s_m=0.0_wp !For H2O, sm is zero
      end if
      
   end subroutine init_smd_h2o

   subroutine init_smd_ot(param,n,alpha,beta,msurft,arom,fclbr)
      !> Refraction index, Abrahams hydrogen bond accidity and Abrahams hydrogen bond basicity of the Solvent
      real(wp), intent(in) :: n, alpha, beta
      !> macroscopic surface tension, fraction of aromatic atoms and fraction of f,cl and br of the solvent
      real(wp), intent(in) :: msurft,arom,fclbr
      !> smd Parameters
      type(smd_param), intent(out) :: param

      !> Local Variables
      !> Is there a self defined Parameter file?
      logical :: ex
      !> Laufvariable
      integer :: Z,Z2

      Inquire(file='smd_ot.param',exist=ex)
      
      if (ex) then
         write(*,*) 'Self defined smd Parameters not implemented yet.'
         stop
      else
         Call init_default(.false.)
         do Z=1,max_elem
            param%zk(Z)=ref_zk(1,Z)*n+ref_zk(2,Z)*alpha+ref_zk(3,Z)*beta
            do Z2=1,max_elem
               param%zkk(Z,Z2)=ref_zkk(1,Z,Z2)*n+ref_zkk(2,Z,Z2)*alpha+ref_zkk(3,Z,Z2)*beta
            end do
         end do
         param%s_m=ref_sg*msurft+ref_sr2*(arom**2)+ref_sp2*(fclbr**2)+ref_sb2*(beta**2)
         param%rzkk=ref_rzkk
         param%drzkk=ref_drzkk
         param%nc3=ref_nc3
         param%rnc3=ref_rnc3
         param%drnc3=ref_drnc3
      end if


   end subroutine init_smd_ot

   subroutine init_default(h2o)
      !>Init H2O defaults?
      logical, intent(in) :: h2o

      !Default smd Parameters taken from:
      !Aleksandr V. Marenich, Christopher J. Cramer, and Donald G. Truhlar
      !The Journal of Physical Chemistry B 2009 113 (18), 6378-6396
      !DOI: 10.1021/jp810292n 

      if (h2o) then
         !> Surface Tension for single atoms
         ref_zk_h2o=0.0_wp !Parameters that are not defined are zero
         ref_zk_h2o(1)=48.69_wp !H
         ref_zk_h2o(6)=129.74_wp !C
         ref_zk_h2o(9)=38.18_wp !F
         ref_zk_h2o(17)=9.82_wp !Cl
         ref_zk_h2o(35)=-8.72_wp !Br
         ref_zk_h2o(16)=-9.10 !S
         !> Surface Tension for Atom-Atom Interactions
         ref_zkk_h2o=0.0_wp !Parameters that are not defined are zero
         ref_zkk_h2o(1,6)=-60.77_wp !H,C
         ref_zkk_h2o(6,6)=-72.95_wp !C,C
         ref_zkk_h2o(8,6)=68.69_wp !O,C
         ref_zkk_h2o(7,6)=-48.22_wp !N,C
         ref_zkk_h2o(8,7)=121.98 !O,N
         ref_zkk_h2o(8,15)=68.85 !O,P
         ref_nc3=84.10_wp !Special Parameter NC3
      else
         ref_zk=0.0_wp 
         ref_zkk=0.0_wp !Parameters that are not defined are zero
         !> Refractive Index Parameters
         !> Element dependent
         ref_zk(1,6)=58.10_wp !C
         ref_zk(1,8)=-17.56_wp !O
         ref_zk(1,7)=32.62_wp !N
         ref_zk(1,17)=-24.31_wp !Cl
         ref_zk(1,35)=-35.42_wp !Br
         ref_zk(1,16)=-33.17_wp !S
         ref_zk(1,14)=-18.04_wp !Si
         !> Pair dependent
         ref_zkk(1,1,6)=-36.37_wp !H,C
         ref_zkk(1,6,6)=-62.05_wp !C,C
         ref_zkk(1,1,8)=-19.39_wp !H,O
         ref_zkk(1,8,6)=-15.70_wp !O,C
         ref_zkk(1,6,7)=-99.76_wp !C,N
         !> Abrahams HB acidity parameters
         !> Element dependent
         ref_zk(2,6)=48.10_wp !C
         ref_zk(2,8)=193.06_wp !O
         !> Pair dependent
         ref_zkk(2,8,6)=95.99_wp !O,C
         ref_zkk(2,6,7)=152.20_wp !C,N
         ref_zkk(2,7,6)=-41.00_wp !N,C
         !> Abrahams HB basicity parameters
         !> Element dependent
         ref_zk(3,6)=32.87_wp !C
         ref_zk(3,8)=-43.79_wp !O
         !> Pair dependent
         ref_zkk(3,8,8)=-128.16_wp !O,O
         ref_zkk(3,8,7)=79.13_wp !O,N

         ref_nc3=0.0_wp
      end if

      ref_rzkk=0.0_wp
      ref_drzkk=0.0_wp

      ref_rzkk(1,6)=1.55_wp !H,C
      ref_rzkk(1,8)=1.55_wp !H,O
      ref_rzkk(6,1)=1.55_wp !C,H
      ref_rzkk(6,6)=1.84_wp !C,C
      ref_rzkk(6,7)=1.84_wp !C,N
      ref_rzkk(6,8)=1.84_wp !C,P
      ref_rzkk(6,9)=1.84_wp !C,F
      ref_rzkk(6,15)=2.2_wp !C,P
      ref_rzkk(6,16)=2.2_wp !C,S
      ref_rzkk(6,17)=2.1_wp !C,Cl
      ref_rzkk(6,35)=2.3_wp !C,Br
      ref_rzkk(6,53)=2.6_wp !C,I
      ref_rzkk(7,6)=1.84_wp !N,C
      ref_rzkk(8,6)=1.33_wp !O,C
      ref_rzkk(8,7)=1.5_wp !O,N
      ref_rzkk(8,8)=1.8_wp !O,O
      ref_rzkk(8,15)=2.1_wp !O,P
      ref_rnc3=1.225_wp !Special nc3 Parameter

      ref_drzkk(1,6)=0.3_wp !H,C
      ref_drzkk(1,8)=0.3_wp !H,O
      ref_drzkk(6,1)=0.3_wp !C,H
      ref_drzkk(6,6)=0.3_wp !C,C
      ref_drzkk(6,7)=0.3_wp !C,N
      ref_drzkk(6,8)=0.3_wp !C,P
      ref_drzkk(6,9)=0.3_wp !C,F
      ref_drzkk(6,15)=0.3_wp !C,P
      ref_drzkk(6,16)=0.3_wp !C,S
      ref_drzkk(6,17)=0.3_wp !C,Cl
      ref_drzkk(6,35)=0.3_wp !C,Br
      ref_drzkk(6,53)=0.3_wp !C,I
      ref_drzkk(7,6)=0.3_wp !N,C
      ref_drzkk(8,6)=0.1_wp !O,C
      ref_drzkk(8,7)=0.3_wp !O,N
      ref_drzkk(8,8)=0.3_wp !O,O
      ref_drzkk(8,15)=0.3_wp !O,P
      ref_drnc3=0.065_wp !Special nc3 Parameter
   end subroutine init_default

end module smd_init
