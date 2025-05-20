!------------------------------------------------------------------------------
! © 2021-2023. Triad National Security, LLC. All rights reserved.  This
! program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by Triad
! National Security, LLC for the U.S.  Department of Energy/National
! Nuclear Security Administration. All rights in the program are
! reserved by Triad National Security, LLC, and the U.S. Department of
! Energy/National Nuclear Security Administration. The Government is
! granted for itself and others acting on its behalf a nonexclusive,
! paid-up, irrevocable worldwide license in this material to reproduce,
! prepare derivative works, distribute copies to the public, perform
! publicly and display publicly, and to permit others to do so.
!------------------------------------------------------------------------------

program test_sg_fortran_interface
! modules
use singularity_eos
! no implicit vars
implicit none
! variable declaration
type(sg_eos) :: eos
real(kind=8), dimension(2) :: temps, rhos, sies
logical, dimension(2) :: mask

temps(1) = 1000.d0
temps(2) = 2000.d0

rhos(1) = 0.1d0
rhos(2) = 0.2d0

mask(1) = .false.
mask(2) = .true.

eos = init_IdealGas(1.4d0, 1.0d0)

sies(:) = 0.0d0
sies = eos%InternalEnergyFromDensityTemperature(0.1d0, 3000.d0)
! assert (sies == 3000)

sies(:) = 0.0d0
sies = eos%InternalEnergyFromDensityTemperature(rhos, temps)
! assert sies == (1000, 2000)

sies(:) = 0.0d0
where (mask)
    sies = eos%InternalEnergyFromDensityTemperature(rhos, temps)
endwhere
! assert sies == (0, 2000)

call eos%Finalize()
end program test_sg_fortran_interface
