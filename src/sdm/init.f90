module sdm_init

!--------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------!
!-----               This Module initiates the SDM Parameters               -----!
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
   implicit none
   private
   public :: init_sdm, sdm_param

   integer, parameter :: max_elem = 94
   
   type :: sdm_param
      !> Element dependent Parameters
      real(wp) :: zk(3,max_elem)
      !> Pair/Solvent dependent Parameters
      real(wp) :: zkk(3,max_elem,max_elem)
      !> Pair dependent Parameters
      real(wp) :: rzkk(max_elem,max_elem)
      real(wp) :: drzkk(max_elem,max_elem)
      !> Additional nc3 Parameters
      real(wp) :: nc3
      real(wp) :: rnc3
      real(wp) :: drnc3
   end type sdm_param

 
   real(wp) :: ref_zk_h2o(max_elem)
   
   real(wp) :: ref_zkk_h2o(max_elem,max_elem)

   real(wp) :: ref_zk(3,max_elem)

   real(wp) :: ref_zkk(3,max_elem,max_elem)

   real(wp) :: ref_rzkk(max_elem,max_elem), ref_drzkk(max_elem, max_elem)

   real(wp) :: ref_nc3, ref_rnc3, ref_drnc3

   interface init_sdm
      module procedure :: init_sdm_h2o
   !   module procedure :: init_sdm
   end interface

contains

   subroutine init_sdm_h2o(param,symbol,zk_in,zkk_in,nc3_in)
      !> Self defined Parameter for the SDM
      real(wp), intent(in), optional :: zk_in(:), zkk_in(:,:), nc3_in
      !> Symbols for self defined Parameters
      character(len=*), intent(in), optional :: symbol
      !> SDM Parameters
      type(sdm_param), intent(out) :: param

      !> Local Variables
      !> Working Parameters
      real(wp) :: zk(max_elem), nc3, rnc3, drnc3
      real(wp), dimension(max_elem,max_elem) :: zkk, rzkk, drzkk
      
      Call init_default(.TRUE.)

      if (present(zk_in)) then
         if (present(zkk_in)) then
            if (present(nc3_in)) then
               write(*,*) "Using self defined Zk, Zkk and nc3 Parameters."
            else
               write(*,*) "Self defined Zk and Zkk Parameters given, but no nc3 Parameters, using default nc3 instead."
               nc3=ref_nc3
            end if
         else
            write(*,*) "Self defined Zk Parameters given, but no Zkk Parameters. &
               &Ignoring the Zk Input and using default SDM Parameters instead."
            zk=ref_zk_h2o
            zkk=ref_zkk_h2o
            nc3=ref_nc3
         end if
      else
         write(*,*) "Using default SDM Parameters."
         param%zk(1,:)=ref_zk_h2o
         param%zkk(1,:,:)=ref_zkk_h2o
         param%rzkk=ref_rzkk
         param%drzkk=ref_drzkk
         param%nc3=ref_nc3
         param%rnc3=ref_rnc3
         param%drnc3=ref_drnc3
      end if
      
   end subroutine init_sdm_h2o

   !subroutine init_sdm(xyz,n,alpha,beta,surft,arom,fclbr,sdm,symbol,zk_in,zkk_in,nc3_in)
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

   subroutine init_default(h2o)
      !>Init H2O defaults?
      logical, intent(in) :: h2o

      !Default SDM Parameters taken from:
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
         write(*,*) "Only H2O implemented!"
         stop
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

end module sdm_init
