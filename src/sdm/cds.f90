module sdm_cds
   use sdm_sigma, only: sdm_surft
   use mctc_env, only: wp
   use mctc_io_symbols, only: to_number
   implicit none

   interface calc_cds
      module procedure :: calc_cds_symbol
      module procedure :: calc_cds_number
   end interface

contains

   subroutine calc_cds_symbol(surft,surface,symbols,id,cds)
      !>SDM Surface Tensions per Atom and for the Solvent
      type(sdm_surft) :: surft
      !>SASA Surface per Atom
      real(wp) :: surface(:)
      !>Identifiers per Atom as Symbols
      character(len=*) :: symbols(:)
      !>CDS Part of the Energy
      real(wp) :: cds
      !>ID of Atom
      integer :: id(:)

      !>CDS Part of the Energy per Atom
      real(wp) :: atom_cds(size(surface))
      !> Laufvariable
      integer :: i
      !>Atomic Number of Atom i
      integer :: Z(size(id))

      do i=1,size(id)
         Z(i)=to_number(symbols(id(i)))
      end do

      cds=0.0_wp
      atom_cds=0.0_wp
      do i=1,size(id)
         atom_cds(i)=surft%sk(Z(i))*surface(i)
      end do

      cds=sum(atom_cds)+surft%sm*sum(surface)
   end subroutine calc_cds_symbol
         
   subroutine calc_cds_number(surft,surface,numbers,cds)
      !>SDM Surface Tensions per Atom and for the Solvent
      type(sdm_surft) :: surft
      !>SASA Surface per Atom
      real(wp) :: surface(:)
      !>Identifiers per Atom as Numbers
      integer :: numbers(:)
      !>CDS Part of the Energy
      real(wp) :: cds

      !>CDS Part of the Energy per Atom
      real(wp) :: atom_cds(size(surface))
      !> Laufvariable
      integer :: i

      cds=0.0_wp
      atom_cds=0.0_wp

      do concurrent (i=1:size(surface))
         atom_cds(i)=surft%sk(numbers(i))*surface(i)
      end do

      cds=sum(atom_cds)+surft%sm*sum(surface)
   end subroutine calc_cds_number


end module sdm_cds
