! This file is part of numsa.
! SPDX-Identifier: LGPL-3.0-or-later
!
! numsa is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! numsa is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with numsa.  If not, see <https://www.gnu.org/licenses/>.

!--------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------!
!-----               This Module initiates the smd Parameters               -----!
!-----            Parameter files need to be appropriately named.           -----!
!-----              h2o: smd_h2o           other solvents: smd_ot           -----!
!-----              Own Parameters need to be added as follows:             -----!
!-----                            #Zk                                       -----!
!-----                            Symbol zk                                 -----!
!-----                            Symbol zk                                 -----!
!-----                            Symbol zk                                 -----!
!-----                            #Zk2 (not for h2o)                        -----!
!-----                            Symbol zk                                 -----!
!-----                            Symbol zk                                 -----!
!-----                            Symbol zk                                 -----!
!-----                            #Zk3 (not for h2o)                        -----!
!-----                            Symbol zk                                 -----!
!-----                            Symbol zk                                 -----!
!-----                            Symbol zk                                 -----!
!-----                            #Zkk Matrix                               -----!
!----- Do not add this line!            Symbol Symbol Symbol                -----!
!-----                            Symbol zkk    zkk    zkk                  -----!
!-----                            Symbol zkk    zkk    zkk                  -----!
!-----                            Symbol zkk    zkk    zkk                  -----! 
!-----                            #Zkk Matrix 2 (not for h2o)               -----!
!-----                            Symbol zkk    zkk    zkk                  -----!
!-----                            Symbol zkk    zkk    zkk                  -----!
!-----                            Symbol zkk    zkk    zkk                  -----! 
!-----                            #Zkk Matrix 3 (not for h2o)               -----!
!-----                            Symbol zkk    zkk    zkk                  -----!
!-----                            Symbol zkk    zkk    zkk                  -----!
!-----                            Symbol zkk    zkk    zkk                  -----! 
!-----                            #rzkk Matrix                              -----!
!-----                            Symbol rzkk   rzkk   rzkk                 -----!
!-----                            Symbol rzkk   rzkk   rzkk                 -----!
!-----                            Symbol rzkk   rzkk   rzkk                 -----! 
!-----                            #drzkk Matrix                              -----!
!-----                            Symbol rzkk   rzkk   rzkk                 -----!
!-----                            Symbol rzkk   rzkk   rzkk                 -----!
!-----                            Symbol rzkk   rzkk   rzkk                 -----! 
!-----                            #NC3                                      -----!
!-----                            nc3                                       -----!
!-----                            rnc3                                      -----!
!-----                            drnc3                                     -----!
!-----                            #Solvent (not for h2o)                    -----!
!-----                            s_g                                       -----!
!-----                            s_r2                                      -----!
!-----                            s_p2                                      -----!
!-----                            s_b2                                      -----!
!--------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------!
module smd_init
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

   real(wp) :: ref_sg

   real(wp) :: ref_sr2

   real(wp) :: ref_sp2

   real(wp) :: ref_sb2

   interface init_smd
      module procedure :: init_smd_def
      module procedure :: init_smd_self
   end interface

contains

   subroutine init_smd_def(param,solvent)
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
         case('orcamet')
            Call init_smd_ot(param,1.3442_wp,0.0700_wp,0.3200_wp,41.25_wp,0.0_wp,0.0_wp)
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

   end subroutine init_smd_def

   subroutine init_smd_self(param,solvent,hd)
      !>Solvent Name
      character(len=*), intent(in) :: solvent
      !>Output Parameters
      type(smd_param), intent(out) :: param
      !>Home Directory for self defined Solvent Parameters
      character(len=*), intent(in) :: hd

      !>Local Variables
      !>Are there Solvent Properties for the Solvent?
      logical :: ex
      !> Variables for Solvent Properties
      real(wp) :: n,alpha,beta,msurft,arom,fclbr

      ex=.false.

      select case(solvent)
         case('h2o', 'water')
            Call init_smd_h2o(param,hd)
         case('methanol','ch4')
            Call init_smd_ot(param,1.3314_wp,0.43_wp,0.47_wp,32.3847_wp,0.0_wp,0.0_wp,hd)
         case('dmso','DMSO')
            Call init_smd_ot(param,1.4772_wp,0.0_wp,0.88_wp,62.6680_wp,0.0_wp,0.0_wp,hd)
         case('acetonitrile')
            Call init_smd_ot(param,1.3421_wp,0.07_wp,0.32_wp,28.4000_wp,0.0_wp,0.0_wp,hd)
         case default
            INQUIRE(file=solvent//".prop",exist=ex)
            if (ex) then
               write(*,*) 'Reading self defined solvent properties for ', solvent,'.'
               Call read_smd(solvent//".prop",n,alpha,beta,msurft,arom,fclbr)
               Call init_smd_ot(param,n,alpha,beta,msurft,arom,fclbr,hd)
            else
               write(*,*) 'No Solvent Properties for ', solvent, '.'
               stop
            end if
      end select

   end subroutine init_smd_self

   subroutine init_smd_h2o(param,hd)
      !> smd Parameters
      type(smd_param), intent(out) :: param
      !> Home Directory for self defined SMD Parameters
      character(len=*),intent(in), optional :: hd

      !> Local Variables
      !> Working Parameters
      real(wp) :: zk(max_elem), nc3, rnc3, drnc3
      real(wp), dimension(max_elem,max_elem) :: zkk, rzkk, drzkk

      !> Is there a self defined Parameter file?
      logical :: ex


      Inquire(file="smd_h2o",exist=ex)
      if (ex) then
         ! write(*,*) "Using smd_h2o from working directory."
         Call read_smd("smd_h2o",ref_zk_h2o,ref_zkk_h2o,ref_rzkk,ref_drzkk,&
               &ref_nc3,ref_rnc3,ref_drnc3)
      else
         if (present(hd)) then
            Inquire(file=hd//"smd_h2o",exist=ex)
               if (ex) then
                  ! write(*,*) "SMD Parameter directory: "//hd//"."
                  Call read_smd(hd//"smd_h2o",ref_zk_h2o,ref_zkk_h2o,ref_rzkk,ref_drzkk,&
                  &ref_nc3,ref_rnc3,ref_drnc3)
               else
                  ! write(*,*) "Using default SMD Parameters."
                  Call init_default(.TRUE.)
               end if
         else
            ! write(*,*) "Using default SMD Parameters."
            Call init_default(.TRUE.)
         end if
      end if

      
      param%zk(:)=ref_zk_h2o
      param%zkk(:,:)=ref_zkk_h2o
      param%rzkk=ref_rzkk
      param%drzkk=ref_drzkk
      param%nc3=ref_nc3
      param%rnc3=ref_rnc3
      param%drnc3=ref_drnc3
      param%s_m=0.0_wp !For H2O, sm is zero
   end subroutine init_smd_h2o

   subroutine init_smd_ot(param,n,alpha,beta,msurft,arom,fclbr,hd)
      !> Refraction index, Abrahams hydrogen bond accidity and Abrahams hydrogen bond basicity of the Solvent
      real(wp), intent(in) :: n, alpha, beta
      !> macroscopic surface tension, fraction of aromatic atoms and fraction of f,cl and br of the solvent
      real(wp), intent(in) :: msurft,arom,fclbr
      !> smd Parameters
      type(smd_param), intent(out) :: param
      !> Home directory for self defined SMD Parameters
      character(len=*),intent(in), optional :: hd

      !> Local Variables
      !> Is there a self defined Parameter file?
      logical :: ex
      !> Laufvariable
      integer :: Z,Z2

      Inquire(file="smd_ot",exist=ex)
      if (ex) then
         ! write(*,*) "Using smd_ot from working directory."
            Call read_smd(hd//"smd_ot",ref_zk,ref_zkk,ref_rzkk,ref_drzkk,&
               &ref_nc3,ref_rnc3,ref_drnc3,ref_sg,ref_sr2,ref_sp2,ref_sb2)
      else
         if (present(hd)) then
            Inquire(file=hd//"smd_ot",exist=ex)
               if (ex) then
                  ! write(*,*) "Using smd_ot from home directory: "//hd//"."
                  Call read_smd(hd//"smd_ot",ref_zk,ref_zkk,ref_rzkk,ref_drzkk,&
                     &ref_nc3,ref_rnc3,ref_drnc3,ref_sg,ref_sr2,ref_sp2,ref_sb2)
               else
                  ! write(*,*) "Using default SMD Parameters."
                  Call init_default(.FALSE.)
               end if
         else
            ! write(*,*) "Using default SMD Parameters."
            Call init_default(.FALSE.)
         end if
      end if

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
         ref_zk_h2o(16)=-9.10_wp !S
         !> Surface Tension for Atom-Atom Interactions
         ref_zkk_h2o=0.0_wp !Parameters that are not defined are zero
         ref_zkk_h2o(1,6)=-60.77_wp !H,C
         ref_zkk_h2o(6,6)=-72.95_wp !C,C
         ref_zkk_h2o(8,6)=68.69_wp !O,C
         ref_zkk_h2o(7,6)=-48.22_wp !N,C
         ref_zkk_h2o(8,7)=121.98_wp !O,N
         ref_zkk_h2o(8,15)=68.85_wp !O,P
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

      ref_sg=0.35_wp
      ref_sr2=-4.19_wp
      ref_sp2=-6.68_wp
      ref_sb2=0.00_wp
   end subroutine init_default
end module smd_init
               
