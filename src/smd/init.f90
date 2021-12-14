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
      real(wp) :: alpha
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

      !>Check for self-defined Parameters first
      INQUIRE(file=solvent//".prop",exist=ex)
      if (ex) then
         write(*,*) 'Reading self defined solvent properties for ', solvent,'.'
         Call read_smd(solvent//".prop",n,alpha,beta,msurft,arom,fclbr)
         Call init_smd_ot(param,n,alpha,beta,msurft,arom,fclbr)
      else
      !> Use default Parameters if no self-defined parameters are found
         select case(solvent)
            case('h2o', 'water')
               Call init_smd_h2o(param)
            case ('methanol')
               Call init_smd_ot(param,1.3288_wp,0.43_wp,0.47_wp,31.77_wp,0.0_wp,0.0_wp)
            case ('2methylpyridine')
               Call init_smd_ot(param,1.4957_wp,0.0_wp,0.58_wp,47.5_wp,0.509796_wp,0.0_wp)
            case ('4methyl2pentanone')
               Call init_smd_ot(param,1.3962_wp,0.0_wp,0.51_wp,33.83_wp,0.0_wp,0.0_wp)
            case ('aceticacid')
               Call init_smd_ot(param,1.372_wp,0.61_wp,0.44_wp,39.01_wp,0.0_wp,0.0_wp)
            case ('acetonitrile')
               Call init_smd_ot(param,1.3442_wp,0.07_wp,0.32_wp,41.25_wp,0.0_wp,0.0_wp)
            case ('acetophenone')
               Call init_smd_ot(param,1.5372_wp,0.0_wp,0.48_wp,56.19_wp,0.444889_wp,0.0_wp)
            case ('aniline')
               Call init_smd_ot(param,1.5863_wp,0.26_wp,0.41_wp,60.62_wp,0.734449_wp,0.0_wp)
            case ('anisole')
               Call init_smd_ot(param,1.5174_wp,0.0_wp,0.29_wp,50.52_wp,0.5625_wp,0.0_wp)
            case ('benzene')
               Call init_smd_ot(param,1.5011_wp,0.0_wp,0.14_wp,40.62_wp,1.0_wp,0.0_wp)
            case ('benzylalcohol')
               Call init_smd_ot(param,1.5396_wp,0.33_wp,0.56_wp,52.96_wp,0.5625_wp,0.0_wp)
            case ('bromobenzene')
               Call init_smd_ot(param,1.5597_wp,0.0_wp,0.09_wp,50.72_wp,0.734449_wp,0.020449_wp)
            case ('bromoethane')
               Call init_smd_ot(param,1.4239_wp,0.0_wp,0.12_wp,34.0_wp,0.0_wp,0.110889_wp)
            case ('bromoform')
               Call init_smd_ot(param,1.6005_wp,0.15_wp,0.06_wp,64.58_wp,0.0_wp,0.5625_wp)
            case ('bromooctane')
               Call init_smd_ot(param,1.4524_wp,0.0_wp,0.12_wp,41.28_wp,0.0_wp,0.0121_wp)
            case ('butanol')
               Call init_smd_ot(param,1.3993_wp,0.37_wp,0.48_wp,35.88_wp,0.0_wp,0.0_wp)
            case ('butanone')
               Call init_smd_ot(param,1.3788_wp,0.0_wp,0.51_wp,34.5_wp,0.0_wp,0.0_wp)
            case ('butylacetate')
               Call init_smd_ot(param,1.3941_wp,0.0_wp,0.45_wp,35.81_wp,0.0_wp,0.0_wp)
            case ('butylbenzene')
               Call init_smd_ot(param,1.4898_wp,0.0_wp,0.15_wp,41.33_wp,0.36_wp,0.0_wp)
            case ('carbondisulfide')
               Call init_smd_ot(param,1.6319_wp,0.0_wp,0.07_wp,45.45_wp,0.0_wp,0.0_wp)
            case ('carbontet')
               Call init_smd_ot(param,1.4601_wp,0.0_wp,0.0_wp,38.04_wp,0.0_wp,0.64_wp)
            case ('chlorobenzene')
               Call init_smd_ot(param,1.5241_wp,0.0_wp,0.07_wp,47.48_wp,0.734449_wp,0.020449_wp)
            case ('chloroform')
               Call init_smd_ot(param,1.4459_wp,0.15_wp,0.02_wp,38.39_wp,0.0_wp,0.5625_wp)
            case ('chlorohexane')
               Call init_smd_ot(param,1.4199_wp,0.0_wp,0.1_wp,37.03_wp,0.0_wp,0.020449_wp)
            case ('cyclohexane')
               Call init_smd_ot(param,1.4266_wp,0.0_wp,0.0_wp,35.48_wp,0.0_wp,0.0_wp)
            case ('cyclohexanone')
               Call init_smd_ot(param,1.4507_wp,0.0_wp,0.56_wp,49.76_wp,0.0_wp,0.0_wp)
            case ('decalin')
               Call init_smd_ot(param,1.4528_wp,0.0_wp,0.0_wp,43.82_wp,0.0_wp,0.0_wp)
            case ('decane')
               Call init_smd_ot(param,1.4102_wp,0.0_wp,0.0_wp,33.64_wp,0.0_wp,0.0_wp)
            case ('decanol')
               Call init_smd_ot(param,1.4372_wp,0.37_wp,0.48_wp,41.04_wp,0.0_wp,0.0_wp)
            case ('dibromoethane')
               Call init_smd_ot(param,1.5387_wp,0.1_wp,0.17_wp,56.93_wp,0.0_wp,0.25_wp)
            case ('dibutylether')
               Call init_smd_ot(param,1.3992_wp,0.0_wp,0.45_wp,35.98_wp,0.0_wp,0.0_wp)
            case ('dichloroethane')
               Call init_smd_ot(param,1.4448_wp,0.1_wp,0.11_wp,45.86_wp,0.0_wp,0.25_wp)
            case ('diethylether')
               Call init_smd_ot(param,1.3526_wp,0.0_wp,0.41_wp,23.96_wp,0.0_wp,0.0_wp)
            case ('diisopropylether')
               Call init_smd_ot(param,1.3679_wp,0.0_wp,0.41_wp,24.86_wp,0.0_wp,0.0_wp)
            case ('dimethylacetamide')
               Call init_smd_ot(param,1.438_wp,0.0_wp,0.78_wp,47.62_wp,0.0_wp,0.0_wp)
            case ('dimethylformamide')
               Call init_smd_ot(param,1.4305_wp,0.0_wp,0.74_wp,49.56_wp,0.0_wp,0.0_wp)
            case ('dimethylpyridine')
               Call init_smd_ot(param,1.4953_wp,0.0_wp,0.63_wp,44.64_wp,0.390625_wp,0.0_wp)
            case ('dimethylsulfoxide','dmso','DMSO')
               Call init_smd_ot(param,1.417_wp,0.0_wp,0.88_wp,61.78_wp,0.0_wp,0.0_wp)
            case ('dodecane')
               Call init_smd_ot(param,1.4216_wp,0.0_wp,0.0_wp,35.85_wp,0.0_wp,0.0_wp)
            case ('ethanol')
               Call init_smd_ot(param,1.3611_wp,0.37_wp,0.48_wp,31.62_wp,0.0_wp,0.0_wp)
            case ('ethoxybenzene')
               Call init_smd_ot(param,1.5076_wp,0.0_wp,0.32_wp,46.65_wp,0.444889_wp,0.0_wp)
            case ('ethylacetate')
               Call init_smd_ot(param,1.3723_wp,0.0_wp,0.45_wp,33.67_wp,0.0_wp,0.0_wp)
            case ('ethylbenzene')
               Call init_smd_ot(param,1.4959_wp,0.0_wp,0.15_wp,41.38_wp,0.5625_wp,0.0_wp)
            case ('fluorobenzene')
               Call init_smd_ot(param,1.4684_wp,0.0_wp,0.1_wp,38.37_wp,0.734449_wp,0.020449_wp)
            case ('fluoroctane')
               Call init_smd_ot(param,1.3935_wp,0.0_wp,0.1_wp,33.92_wp,0.0_wp,0.012321_wp)
            case ('heptane')
               Call init_smd_ot(param,1.3878_wp,0.0_wp,0.0_wp,28.28_wp,0.0_wp,0.0_wp)
            case ('hexadecane')
               Call init_smd_ot(param,1.4345_wp,0.0_wp,0.0_wp,38.93_wp,0.0_wp,0.0_wp)
            case ('hexadecyliodide')
               Call init_smd_ot(param,1.4806_wp,0.0_wp,0.15_wp,46.48_wp,0.0_wp,0.0_wp)
            case ('hexane')
               Call init_smd_ot(param,1.3749_wp,0.0_wp,0.0_wp,25.7495_wp,0.0_wp,0.0_wp)
            case ('hexanol')
               Call init_smd_ot(param,1.4178_wp,0.37_wp,0.48_wp,37.15_wp,0.0_wp,0.0_wp)
            case ('iodobenzene')
               Call init_smd_ot(param,1.62_wp,0.0_wp,0.12_wp,55.72_wp,0.734449_wp,0.0_wp)
            case ('isobutanol')
               Call init_smd_ot(param,1.3955_wp,0.37_wp,0.48_wp,32.38_wp,0.0_wp,0.0_wp)
            case ('isooctane')
               Call init_smd_ot(param,1.3915_wp,0.0_wp,0.0_wp,26.38_wp,0.0_wp,0.0_wp)
            case ('isopropanol')
               Call init_smd_ot(param,1.3776_wp,0.33_wp,0.56_wp,30.13_wp,0.0_wp,0.0_wp)
            case ('isopropylbenzene')
               Call init_smd_ot(param,1.4915_wp,0.0_wp,0.16_wp,39.85_wp,0.444889_wp,0.0_wp)
            case ('isopropyltoluene')
               Call init_smd_ot(param,1.4909_wp,0.0_wp,0.19_wp,38.34_wp,0.36_wp,0.0_wp)
            case ('mcresol')
               Call init_smd_ot(param,1.5438_wp,0.57_wp,0.34_wp,51.37_wp,0.5625_wp,0.0_wp)
            case ('mesitylene')
               Call init_smd_ot(param,1.4994_wp,0.0_wp,0.19_wp,39.65_wp,0.444889_wp,0.0_wp)
            case ('methoxyethanol')
               Call init_smd_ot(param,1.4024_wp,0.3_wp,0.84_wp,44.39_wp,0.0_wp,0.0_wp)
            case ('methylenechloride')
               Call init_smd_ot(param,1.4242_wp,0.1_wp,0.05_wp,39.15_wp,0.0_wp,0.444889_wp)
            case ('methylformamide')
               Call init_smd_ot(param,1.4319_wp,0.4_wp,0.55_wp,55.4372_wp,0.0_wp,0.0_wp)
            case ('nitrobenzene')
               Call init_smd_ot(param,1.5562_wp,0.0_wp,0.28_wp,57.54_wp,0.444889_wp,0.0_wp)
            case ('nitroethane')
               Call init_smd_ot(param,1.3917_wp,0.02_wp,0.33_wp,46.25_wp,0.0_wp,0.0_wp)
            case ('nitromethane')
               Call init_smd_ot(param,1.3817_wp,0.06_wp,0.31_wp,52.58_wp,0.0_wp,0.0_wp)
            case ('nonane')
               Call init_smd_ot(param,1.4054_wp,0.0_wp,0.0_wp,32.21_wp,0.0_wp,0.0_wp)
            case ('nonanol')
               Call init_smd_ot(param,1.4333_wp,0.37_wp,0.48_wp,40.14_wp,0.0_wp,0.0_wp)
            case ('octane')
               Call init_smd_ot(param,1.3974_wp,0.0_wp,0.0_wp,30.4273_wp,0.0_wp,0.0_wp)
            case ('octanol')
               Call init_smd_ot(param,1.4295_wp,0.37_wp,0.48_wp,39.01_wp,0.0_wp,0.0_wp)
            case ('odichlorobenzene')
               Call init_smd_ot(param,1.5515_wp,0.0_wp,0.04_wp,52.72_wp,0.5625_wp,0.0625_wp)
            case ('onitrotoluene')
               Call init_smd_ot(param,1.545_wp,0.0_wp,0.27_wp,59.12_wp,0.36_wp,0.0_wp)
            case ('pentadecane')
               Call init_smd_ot(param,1.4315_wp,0.0_wp,0.0_wp,38.34_wp,0.0_wp,0.0_wp)
            case ('pentane')
               Call init_smd_ot(param,1.3575_wp,0.0_wp,0.0_wp,22.3_wp,0.0_wp,0.0_wp)
            case ('pentanol')
               Call init_smd_ot(param,1.4101_wp,0.37_wp,0.48_wp,36.5_wp,0.0_wp,0.0_wp)
            case ('perfluorobenzene')
               Call init_smd_ot(param,1.3777_wp,0.0_wp,0.0_wp,31.74_wp,0.25_wp,0.25_wp)
            case ('phenylether')
               Call init_smd_ot(param,1.5787_wp,0.0_wp,0.2_wp,38.5_wp,0.851929_wp,0.0_wp)
            case ('propanol')
               Call init_smd_ot(param,1.385_wp,0.37_wp,0.48_wp,33.57_wp,0.0_wp,0.0_wp)
            case ('pyridine')
               Call init_smd_ot(param,1.5095_wp,0.0_wp,0.52_wp,52.62_wp,0.693889_wp,0.0_wp)
            case ('secbutanol')
               Call init_smd_ot(param,1.3978_wp,0.33_wp,0.56_wp,32.44_wp,0.0_wp,0.0_wp)
            case ('secbutylbenzene')
               Call init_smd_ot(param,1.4895_wp,0.0_wp,0.16_wp,40.35_wp,0.36_wp,0.0_wp)
            case ('tbutylbenzene')
               Call init_smd_ot(param,1.4927_wp,0.0_wp,0.16_wp,39.78_wp,0.36_wp,0.0_wp)
            case ('tetrachloroethene')
               Call init_smd_ot(param,1.5053_wp,0.0_wp,0.0_wp,45.19_wp,0.0_wp,0.444889_wp)
            case ('tetrahydrofuran')
               Call init_smd_ot(param,1.405_wp,0.0_wp,0.48_wp,39.44_wp,0.0_wp,0.0_wp)
            case ('tetrahydrothiophenedioxide')
               Call init_smd_ot(param,1.4833_wp,0.0_wp,0.88_wp,87.49_wp,0.0_wp,0.0_wp)
            case ('tetralin')
               Call init_smd_ot(param,1.5413_wp,0.0_wp,0.19_wp,47.74_wp,0.36_wp,0.0_wp)
            case ('toluene')
               Call init_smd_ot(param,1.4961_wp,0.0_wp,0.14_wp,40.2_wp,0.734449_wp,0.0_wp)
            case ('tributylphosphate')
               Call init_smd_ot(param,1.4224_wp,0.0_wp,1.21_wp,27.55_wp,0.0_wp,0.0_wp)
            case ('triethylamine')
               Call init_smd_ot(param,1.401_wp,0.0_wp,0.79_wp,29.1_wp,0.0_wp,0.0_wp)
            case ('trimethylbenzene')
               Call init_smd_ot(param,1.5048_wp,0.0_wp,0.19_wp,42.03_wp,0.444889_wp,0.0_wp)
            case ('undecane')
               Call init_smd_ot(param,1.4398_wp,0.0_wp,0.0_wp,34.85_wp,0.0_wp,0.0_wp)
            case ('xylene')
               Call init_smd_ot(param,1.4995_wp,0.0_wp,0.16_wp,41.38_wp,0.5625_wp,0.0_wp)
            case default
               write(*,*) 'No Solvent Properties for ', solvent, '.'
               stop
         end select
      end if
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

      !> Check for self-defined parameters first.
      INQUIRE(file=solvent//".prop",exist=ex)
      if (ex) then
         write(*,*) 'Reading self defined solvent properties for ', solvent,'.'
         Call read_smd(solvent//".prop",n,alpha,beta,msurft,arom,fclbr)
         Call init_smd_ot(param,n,alpha,beta,msurft,arom,fclbr,trim(hd))
      else
      !> Use default parameters if no self-defined Parameters are found
         select case(solvent)
            case('h2o', 'water')
               Call init_smd_h2o(param,trim(hd))
            case ('methanol')
               Call init_smd_ot(param,1.3288_wp,0.43_wp,0.47_wp,31.77_wp,0.0_wp,0.0_wp,trim(hd))
            case ('2methylpyridine')
               Call init_smd_ot(param,1.4957_wp,0.0_wp,0.58_wp,47.5_wp,0.509796_wp,0.0_wp,trim(hd))
            case ('4methyl2pentanone')
               Call init_smd_ot(param,1.3962_wp,0.0_wp,0.51_wp,33.83_wp,0.0_wp,0.0_wp,trim(hd))
            case ('aceticacid')
               Call init_smd_ot(param,1.372_wp,0.61_wp,0.44_wp,39.01_wp,0.0_wp,0.0_wp,trim(hd))
            case ('acetonitrile')
               Call init_smd_ot(param,1.3442_wp,0.07_wp,0.32_wp,41.25_wp,0.0_wp,0.0_wp,trim(hd))
            case ('acetophenone')
               Call init_smd_ot(param,1.5372_wp,0.0_wp,0.48_wp,56.19_wp,0.444889_wp,0.0_wp,trim(hd))
            case ('aniline')
               Call init_smd_ot(param,1.5863_wp,0.26_wp,0.41_wp,60.62_wp,0.734449_wp,0.0_wp,trim(hd))
            case ('anisole')
               Call init_smd_ot(param,1.5174_wp,0.0_wp,0.29_wp,50.52_wp,0.5625_wp,0.0_wp,trim(hd))
            case ('benzene')
               Call init_smd_ot(param,1.5011_wp,0.0_wp,0.14_wp,40.62_wp,1.0_wp,0.0_wp,trim(hd))
            case ('benzylalcohol')
               Call init_smd_ot(param,1.5396_wp,0.33_wp,0.56_wp,52.96_wp,0.5625_wp,0.0_wp,trim(hd))
            case ('bromobenzene')
               Call init_smd_ot(param,1.5597_wp,0.0_wp,0.09_wp,50.72_wp,0.734449_wp,0.020449_wp,trim(hd))
            case ('bromoethane')
               Call init_smd_ot(param,1.4239_wp,0.0_wp,0.12_wp,34.0_wp,0.0_wp,0.110889_wp,trim(hd))
            case ('bromoform')
               Call init_smd_ot(param,1.6005_wp,0.15_wp,0.06_wp,64.58_wp,0.0_wp,0.5625_wp,trim(hd))
            case ('bromooctane')
               Call init_smd_ot(param,1.4524_wp,0.0_wp,0.12_wp,41.28_wp,0.0_wp,0.0121_wp,trim(hd))
            case ('butanol')
               Call init_smd_ot(param,1.3993_wp,0.37_wp,0.48_wp,35.88_wp,0.0_wp,0.0_wp,trim(hd))
            case ('butanone')
               Call init_smd_ot(param,1.3788_wp,0.0_wp,0.51_wp,34.5_wp,0.0_wp,0.0_wp,trim(hd))
            case ('butylacetate')
               Call init_smd_ot(param,1.3941_wp,0.0_wp,0.45_wp,35.81_wp,0.0_wp,0.0_wp,trim(hd))
            case ('butylbenzene')
               Call init_smd_ot(param,1.4898_wp,0.0_wp,0.15_wp,41.33_wp,0.36_wp,0.0_wp,trim(hd))
            case ('carbondisulfide')
               Call init_smd_ot(param,1.6319_wp,0.0_wp,0.07_wp,45.45_wp,0.0_wp,0.0_wp,trim(hd))
            case ('carbontet')
               Call init_smd_ot(param,1.4601_wp,0.0_wp,0.0_wp,38.04_wp,0.0_wp,0.64_wp,trim(hd))
            case ('chlorobenzene')
               Call init_smd_ot(param,1.5241_wp,0.0_wp,0.07_wp,47.48_wp,0.734449_wp,0.020449_wp,trim(hd))
            case ('chloroform')
               Call init_smd_ot(param,1.4459_wp,0.15_wp,0.02_wp,38.39_wp,0.0_wp,0.5625_wp,trim(hd))
            case ('chlorohexane')
               Call init_smd_ot(param,1.4199_wp,0.0_wp,0.1_wp,37.03_wp,0.0_wp,0.020449_wp,trim(hd))
            case ('cyclohexane')
               Call init_smd_ot(param,1.4266_wp,0.0_wp,0.0_wp,35.48_wp,0.0_wp,0.0_wp,trim(hd))
            case ('cyclohexanone')
               Call init_smd_ot(param,1.4507_wp,0.0_wp,0.56_wp,49.76_wp,0.0_wp,0.0_wp,trim(hd))
            case ('decalin')
               Call init_smd_ot(param,1.4528_wp,0.0_wp,0.0_wp,43.82_wp,0.0_wp,0.0_wp,trim(hd))
            case ('decane')
               Call init_smd_ot(param,1.4102_wp,0.0_wp,0.0_wp,33.64_wp,0.0_wp,0.0_wp,trim(hd))
            case ('decanol')
               Call init_smd_ot(param,1.4372_wp,0.37_wp,0.48_wp,41.04_wp,0.0_wp,0.0_wp,trim(hd))
            case ('dibromoethane')
               Call init_smd_ot(param,1.5387_wp,0.1_wp,0.17_wp,56.93_wp,0.0_wp,0.25_wp,trim(hd))
            case ('dibutylether')
               Call init_smd_ot(param,1.3992_wp,0.0_wp,0.45_wp,35.98_wp,0.0_wp,0.0_wp,trim(hd))
            case ('dichloroethane')
               Call init_smd_ot(param,1.4448_wp,0.1_wp,0.11_wp,45.86_wp,0.0_wp,0.25_wp,trim(hd))
            case ('diethylether')
               Call init_smd_ot(param,1.3526_wp,0.0_wp,0.41_wp,23.96_wp,0.0_wp,0.0_wp,trim(hd))
            case ('diisopropylether')
               Call init_smd_ot(param,1.3679_wp,0.0_wp,0.41_wp,24.86_wp,0.0_wp,0.0_wp,trim(hd))
            case ('dimethylacetamide')
               Call init_smd_ot(param,1.438_wp,0.0_wp,0.78_wp,47.62_wp,0.0_wp,0.0_wp,trim(hd))
            case ('dimethylformamide')
               Call init_smd_ot(param,1.4305_wp,0.0_wp,0.74_wp,49.56_wp,0.0_wp,0.0_wp,trim(hd))
            case ('dimethylpyridine')
               Call init_smd_ot(param,1.4953_wp,0.0_wp,0.63_wp,44.64_wp,0.390625_wp,0.0_wp,trim(hd))
            case ('dimethylsulfoxide','dmso','DMSO')
               Call init_smd_ot(param,1.417_wp,0.0_wp,0.88_wp,61.78_wp,0.0_wp,0.0_wp,trim(hd))
            case ('dodecane')
               Call init_smd_ot(param,1.4216_wp,0.0_wp,0.0_wp,35.85_wp,0.0_wp,0.0_wp,trim(hd))
            case ('ethanol')
               Call init_smd_ot(param,1.3611_wp,0.37_wp,0.48_wp,31.62_wp,0.0_wp,0.0_wp,trim(hd))
            case ('ethoxybenzene')
               Call init_smd_ot(param,1.5076_wp,0.0_wp,0.32_wp,46.65_wp,0.444889_wp,0.0_wp,trim(hd))
            case ('ethylacetate')
               Call init_smd_ot(param,1.3723_wp,0.0_wp,0.45_wp,33.67_wp,0.0_wp,0.0_wp,trim(hd))
            case ('ethylbenzene')
               Call init_smd_ot(param,1.4959_wp,0.0_wp,0.15_wp,41.38_wp,0.5625_wp,0.0_wp,trim(hd))
            case ('fluorobenzene')
               Call init_smd_ot(param,1.4684_wp,0.0_wp,0.1_wp,38.37_wp,0.734449_wp,0.020449_wp,trim(hd))
            case ('fluoroctane')
               Call init_smd_ot(param,1.3935_wp,0.0_wp,0.1_wp,33.92_wp,0.0_wp,0.012321_wp,trim(hd))
            case ('heptane')
               Call init_smd_ot(param,1.3878_wp,0.0_wp,0.0_wp,28.28_wp,0.0_wp,0.0_wp,trim(hd))
            case ('hexadecane')
               Call init_smd_ot(param,1.4345_wp,0.0_wp,0.0_wp,38.93_wp,0.0_wp,0.0_wp,trim(hd))
            case ('hexadecyliodide')
               Call init_smd_ot(param,1.4806_wp,0.0_wp,0.15_wp,46.48_wp,0.0_wp,0.0_wp,trim(hd))
            case ('hexane')
               Call init_smd_ot(param,1.3749_wp,0.0_wp,0.0_wp,25.7495_wp,0.0_wp,0.0_wp,trim(hd))
            case ('hexanol')
               Call init_smd_ot(param,1.4178_wp,0.37_wp,0.48_wp,37.15_wp,0.0_wp,0.0_wp,trim(hd))
            case ('iodobenzene')
               Call init_smd_ot(param,1.62_wp,0.0_wp,0.12_wp,55.72_wp,0.734449_wp,0.0_wp,trim(hd))
            case ('isobutanol')
               Call init_smd_ot(param,1.3955_wp,0.37_wp,0.48_wp,32.38_wp,0.0_wp,0.0_wp,trim(hd))
            case ('isooctane')
               Call init_smd_ot(param,1.3915_wp,0.0_wp,0.0_wp,26.38_wp,0.0_wp,0.0_wp,trim(hd))
            case ('isopropanol')
               Call init_smd_ot(param,1.3776_wp,0.33_wp,0.56_wp,30.13_wp,0.0_wp,0.0_wp,trim(hd))
            case ('isopropylbenzene')
               Call init_smd_ot(param,1.4915_wp,0.0_wp,0.16_wp,39.85_wp,0.444889_wp,0.0_wp,trim(hd))
            case ('isopropyltoluene')
               Call init_smd_ot(param,1.4909_wp,0.0_wp,0.19_wp,38.34_wp,0.36_wp,0.0_wp,trim(hd))
            case ('mcresol')
               Call init_smd_ot(param,1.5438_wp,0.57_wp,0.34_wp,51.37_wp,0.5625_wp,0.0_wp,trim(hd))
            case ('mesitylene')
               Call init_smd_ot(param,1.4994_wp,0.0_wp,0.19_wp,39.65_wp,0.444889_wp,0.0_wp,trim(hd))
            case ('methoxyethanol')
               Call init_smd_ot(param,1.4024_wp,0.3_wp,0.84_wp,44.39_wp,0.0_wp,0.0_wp,trim(hd))
            case ('methylenechloride')
               Call init_smd_ot(param,1.4242_wp,0.1_wp,0.05_wp,39.15_wp,0.0_wp,0.444889_wp,trim(hd))
            case ('methylformamide')
               Call init_smd_ot(param,1.4319_wp,0.4_wp,0.55_wp,55.4372_wp,0.0_wp,0.0_wp,trim(hd))
            case ('nitrobenzene')
               Call init_smd_ot(param,1.5562_wp,0.0_wp,0.28_wp,57.54_wp,0.444889_wp,0.0_wp,trim(hd))
            case ('nitroethane')
               Call init_smd_ot(param,1.3917_wp,0.02_wp,0.33_wp,46.25_wp,0.0_wp,0.0_wp,trim(hd))
            case ('nitromethane')
               Call init_smd_ot(param,1.3817_wp,0.06_wp,0.31_wp,52.58_wp,0.0_wp,0.0_wp,trim(hd))
            case ('nonane')
               Call init_smd_ot(param,1.4054_wp,0.0_wp,0.0_wp,32.21_wp,0.0_wp,0.0_wp,trim(hd))
            case ('nonanol')
               Call init_smd_ot(param,1.4333_wp,0.37_wp,0.48_wp,40.14_wp,0.0_wp,0.0_wp,trim(hd))
            case ('octane')
               Call init_smd_ot(param,1.3974_wp,0.0_wp,0.0_wp,30.4273_wp,0.0_wp,0.0_wp,trim(hd))
            case ('octanol')
               Call init_smd_ot(param,1.4295_wp,0.37_wp,0.48_wp,39.01_wp,0.0_wp,0.0_wp,trim(hd))
            case ('odichlorobenzene')
               Call init_smd_ot(param,1.5515_wp,0.0_wp,0.04_wp,52.72_wp,0.5625_wp,0.0625_wp,trim(hd))
            case ('onitrotoluene')
               Call init_smd_ot(param,1.545_wp,0.0_wp,0.27_wp,59.12_wp,0.36_wp,0.0_wp,trim(hd))
            case ('pentadecane')
               Call init_smd_ot(param,1.4315_wp,0.0_wp,0.0_wp,38.34_wp,0.0_wp,0.0_wp,trim(hd))
            case ('pentane')
               Call init_smd_ot(param,1.3575_wp,0.0_wp,0.0_wp,22.3_wp,0.0_wp,0.0_wp,trim(hd))
            case ('pentanol')
               Call init_smd_ot(param,1.4101_wp,0.37_wp,0.48_wp,36.5_wp,0.0_wp,0.0_wp,trim(hd))
            case ('perfluorobenzene')
               Call init_smd_ot(param,1.3777_wp,0.0_wp,0.0_wp,31.74_wp,0.25_wp,0.25_wp,trim(hd))
            case ('phenylether')
               Call init_smd_ot(param,1.5787_wp,0.0_wp,0.2_wp,38.5_wp,0.851929_wp,0.0_wp,trim(hd))
            case ('propanol')
               Call init_smd_ot(param,1.385_wp,0.37_wp,0.48_wp,33.57_wp,0.0_wp,0.0_wp,trim(hd))
            case ('pyridine')
               Call init_smd_ot(param,1.5095_wp,0.0_wp,0.52_wp,52.62_wp,0.693889_wp,0.0_wp,trim(hd))
            case ('secbutanol')
               Call init_smd_ot(param,1.3978_wp,0.33_wp,0.56_wp,32.44_wp,0.0_wp,0.0_wp,trim(hd))
            case ('secbutylbenzene')
               Call init_smd_ot(param,1.4895_wp,0.0_wp,0.16_wp,40.35_wp,0.36_wp,0.0_wp,trim(hd))
            case ('tbutylbenzene')
               Call init_smd_ot(param,1.4927_wp,0.0_wp,0.16_wp,39.78_wp,0.36_wp,0.0_wp,trim(hd))
            case ('tetrachloroethene')
               Call init_smd_ot(param,1.5053_wp,0.0_wp,0.0_wp,45.19_wp,0.0_wp,0.444889_wp,trim(hd))
            case ('tetrahydrofuran')
               Call init_smd_ot(param,1.405_wp,0.0_wp,0.48_wp,39.44_wp,0.0_wp,0.0_wp,trim(hd))
            case ('tetrahydrothiophenedioxide')
               Call init_smd_ot(param,1.4833_wp,0.0_wp,0.88_wp,87.49_wp,0.0_wp,0.0_wp,trim(hd))
            case ('tetralin')
               Call init_smd_ot(param,1.5413_wp,0.0_wp,0.19_wp,47.74_wp,0.36_wp,0.0_wp,trim(hd))
            case ('toluene')
               Call init_smd_ot(param,1.4961_wp,0.0_wp,0.14_wp,40.2_wp,0.734449_wp,0.0_wp,trim(hd))
            case ('tributylphosphate')
               Call init_smd_ot(param,1.4224_wp,0.0_wp,1.21_wp,27.55_wp,0.0_wp,0.0_wp,trim(hd))
            case ('triethylamine')
               Call init_smd_ot(param,1.401_wp,0.0_wp,0.79_wp,29.1_wp,0.0_wp,0.0_wp,trim(hd))
            case ('trimethylbenzene')
               Call init_smd_ot(param,1.5048_wp,0.0_wp,0.19_wp,42.03_wp,0.444889_wp,0.0_wp,trim(hd))
            case ('undecane')
               Call init_smd_ot(param,1.4398_wp,0.0_wp,0.0_wp,34.85_wp,0.0_wp,0.0_wp,trim(hd))
            case ('xylene')
               Call init_smd_ot(param,1.4995_wp,0.0_wp,0.16_wp,41.38_wp,0.5625_wp,0.0_wp,trim(hd))
            case default
               write(*,*) 'No Solvent Properties for ', solvent, '.'
               stop
         end select
      end if

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
                  INQUIRE(file=hd,exist=ex)
                  if (ex) then
                     Call read_smd(hd,ref_zk_h2o,ref_zkk_h2o,ref_rzkk,ref_drzkk,&
                     &ref_nc3,ref_rnc3,ref_drnc3)
                  else
                     ! write(*,*) "Using default SMD Parameters."
                     Call init_default(.TRUE.)
                  end if
               end if
         else
            ! write(*,*) "Using default SMD Parameters."
            Call init_default(.TRUE.)
         end if
      end if

      param%alpha=0.82_wp ! Alpha is needed for O-Radius scaling 
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
                  INQUIRE(file=hd,exist=ex)
                  if (ex) then
                     Call read_smd(hd,ref_zk,ref_zkk,ref_rzkk,ref_drzkk,&
                     &ref_nc3,ref_rnc3,ref_drnc3,ref_sg,ref_sr2,ref_sp2,ref_sb2)
                  else
                  ! write(*,*) "Using default SMD Parameters."
                     Call init_default(.FALSE.)
                  end if
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
      param%alpha=alpha ! Alpha is needed for O Radius scaling

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
               
