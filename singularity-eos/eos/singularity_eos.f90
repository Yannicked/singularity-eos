module singularity_eos
   use, intrinsic :: iso_c_binding
   implicit none

   public :: init_IdealGas, init_SpinerDependsRhoT, sg_eos

   type :: sg_eos
      ! Pointer to C++ EOS variant
      type(c_ptr), private :: ptr = C_NULL_PTR

   contains
      ! EOS API
      procedure :: InternalEnergyFromDensityTemperature
      procedure :: PressureFromDensityTemperature
      procedure :: BulkModulusFromDensityTemperature
      procedure :: SpecificHeatFromDensityTemperature

      ! EOS Properties
      procedure :: MeanAtomicMass
      procedure :: MeanAtomicNumber

      procedure :: Finalize
   end type sg_eos

   interface
      function init_sg_IdealGas(gamma_minus_1, c_v) &
         bind(C, name='init_sg_IdealGas') result(eos)
         import
         real(kind=c_double), value, intent(in) :: gamma_minus_1, c_v
         type(c_ptr) :: eos
      end function init_sg_IdealGas
   end interface

   interface
      function init_sg_SpinerDependsRhoT(filename, id, split) &
         bind(C, name='init_sg_SpinerDependsRhoT') result(eos)
         import
         character(kind=c_char), intent(in) :: filename(*)
         integer(c_int), value, intent(in) :: id, split
         type(c_ptr) :: eos
      end function init_sg_SpinerDependsRhoT
   end interface

   interface
      pure function get_sg_InternalEnergyFromDensityTemperature(eos, rho, temperature) &
         bind(C, name='get_sg_InternalEnergyFromDensityTemperature') result(sie)
         import
         type(c_ptr), value, intent(in) :: eos
         real(kind=c_double), value, intent(in) :: rho, temperature
         real(kind=c_double) :: sie
      end function
   end interface

   interface
      pure function get_sg_PressureFromDensityTemperature(eos, rho, temperature) &
         bind(C, name='get_sg_PressureFromDensityTemperature') result(pressure)
         import
         type(c_ptr), value, intent(in) :: eos
         real(kind=c_double), value, intent(in) :: rho, temperature
         real(kind=c_double) :: pressure
      end function
   end interface

   interface
      pure function get_sg_BulkModulusFromDensityTemperature(eos, rho, temperature) &
         bind(C, name='get_sg_BulkModulusFromDensityTemperature') result(bulk_modulus)
         import
         type(c_ptr), value, intent(in) :: eos
         real(kind=c_double), value, intent(in) :: rho, temperature
         real(kind=c_double) :: bulk_modulus
      end function
   end interface

   interface
      pure function get_sg_SpecificHeatFromDensityTemperature(eos, rho, temperature) &
         bind(C, name='get_sg_SpecificHeatFromDensityTemperature') result(specific_heat)
         import
         type(c_ptr), value, intent(in) :: eos
         real(kind=c_double), value, intent(in) :: rho, temperature
         real(kind=c_double) :: specific_heat
      end function
   end interface

   interface
      pure function get_sg_MeanAtomicMass(eos) &
         bind(C, name='get_sg_MeanAtomicMass') result(a_bar)
         import
         type(c_ptr), value, intent(in) :: eos
         real(kind=c_double) :: a_bar
      end function
   end interface

   interface
      pure function get_sg_MeanAtomicNumber(eos) &
         bind(C, name='get_sg_MeanAtomicNumber') result(z_bar)
         import
         type(c_ptr), value, intent(in) :: eos
         real(kind=c_double) :: z_bar
      end function
   end interface

   interface
      subroutine finalize_sg_eos(eos) bind(C, name='finalize_sg_eos')
         import
         type(c_ptr), value, intent(in) :: eos
      end subroutine
   end interface
contains
   function init_IdealGas(gamma_minus_1, c_v) result(self)
      real(kind=c_double), intent(in) :: gamma_minus_1, c_v
      type(sg_eos) :: self

      self%ptr = init_sg_IdealGas(gamma_minus_1, c_v)
   end function

   function init_SpinerDependsRhoT(filename, id, split) result(self)
      character(len=*, kind=c_char), intent(in) :: filename
      integer(c_int), intent(in) :: id
      integer(c_int), optional, intent(in) :: split
      type(sg_eos) :: self

      integer(c_int) :: split_use

      if (present(split)) then
         split_use = split
      else
         split_use = 0
      end if

      self%ptr = init_sg_SpinerDependsRhoT(trim(filename)//C_NULL_CHAR, id, split_use)
   end function

   elemental function InternalEnergyFromDensityTemperature(self, rho, temperature) result(sie)
      class(sg_eos), intent(in) :: self
      real(kind=c_double), intent(in) :: rho, temperature
      real(kind=c_double) :: sie

      sie = get_sg_InternalEnergyFromDensityTemperature(self%ptr, rho, temperature)
   end function

   elemental function PressureFromDensityTemperature(self, rho, temperature) result(pressure)
      class(sg_eos), intent(in) :: self
      real(kind=c_double), intent(in) :: rho, temperature
      real(kind=c_double) :: pressure

      pressure = get_sg_PressureFromDensityTemperature(self%ptr, rho, temperature)
   end function

   elemental function BulkModulusFromDensityTemperature(self, rho, temperature) result(bulk_modulus)
      class(sg_eos), intent(in) :: self
      real(kind=c_double), intent(in) :: rho, temperature
      real(kind=c_double) :: bulk_modulus

      bulk_modulus = get_sg_BulkModulusFromDensityTemperature(self%ptr, rho, temperature)
   end function

   elemental function SpecificHeatFromDensityTemperature(self, rho, temperature) result(specific_heat)
      class(sg_eos), intent(in) :: self
      real(kind=c_double), intent(in) :: rho, temperature
      real(kind=c_double) :: specific_heat

      specific_heat = get_sg_SpecificHeatFromDensityTemperature(self%ptr, rho, temperature)
   end function

   pure function MeanAtomicMass(self) result(a_bar)
      class(sg_eos), intent(in) :: self
      real(kind=c_double) :: a_bar

      a_bar = get_sg_MeanAtomicMass(self%ptr)
   end function

   pure function MeanAtomicNumber(self) result(z_bar)
      class(sg_eos), intent(in) :: self
      real(kind=c_double) :: z_bar

      z_bar = get_sg_MeanAtomicNumber(self%ptr)
   end function

   subroutine Finalize(self)
      class(sg_eos), intent(in) :: self

      call finalize_sg_eos(self%ptr)
   end subroutine
end module singularity_eos
