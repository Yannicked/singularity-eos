#include <string>

#include <mpark/variant.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_builder.hpp>
#include <singularity-eos/eos/modifiers/ramps_eos.hpp> // Added for BilinearRampEOS
#include <singularity-eos/eos/singularity_eos.hpp>

using namespace singularity;

EOS *init_sg_IdealGas(const double gamma_minus_1, const double c_v) {
  EOS *eos = new EOS;
  *eos = IdealGas(gamma_minus_1, c_v);
  return eos;
}

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5
EOS *init_sg_SpinerDependsRhoT(const char *filename, const int matid,
                               const TableSplit split) {
  EOS *eos = new EOS;
  *eos = SpinerEOSDependsRhoT(std::string(filename), matid, split);
  return eos;
}
#endif

const Real get_sg_InternalEnergyFromDensityTemperature(EOS *eos, const Real density,
                                                       const Real temperature) {
  return eos->InternalEnergyFromDensityTemperature(density, temperature);
}

const Real get_sg_PressureFromDensityTemperature(EOS *eos, const Real density,
                                                 const Real temperature) {
  return eos->PressureFromDensityTemperature(density, temperature);
}

const Real get_sg_BulkModulusFromDensityTemperature(EOS *eos, const Real density,
                                                    const Real temperature) {
  return eos->BulkModulusFromDensityTemperature(density, temperature);
}

const Real get_sg_SpecificHeatFromDensityTemperature(EOS *eos, const Real density,
                                                     const Real temperature) {
  return eos->SpecificHeatFromDensityTemperature(density, temperature);
}

const Real get_sg_MeanAtomicMass(EOS *eos) { return eos->MeanAtomicMass(); }

const Real get_sg_MeanAtomicNumber(EOS *eos) { return eos->MeanAtomicNumber(); }

void modify_sg_ScaledEOS(EOS *eos, const Real scale) {
  *eos = eos->template Modify<ScaledEOS>(scale);
}

void modify_sg_ShiftedEOS(EOS *eos, const Real shift) {
  *eos = eos->template Modify<ShiftedEOS>(shift);
}

void modify_sg_BilinearRampEOS(EOS *eos, const Real r0, const Real a,
                                const Real b, const Real c) {
  // Real r0, a, b, c;
  
  // pAlpha2BilinearRampParams(*eos, alpha0, Pe, Pc, r0, a, b, c);

  *eos = eos->template Modify<BilinearRampEOS>(r0, a, b, c);
}

void finalize_sg_eos(EOS *eos) {
  eos->Finalize();
  delete eos;
  eos = nullptr;
}
