#ifndef SINGULARITY_EOS_HPP
#define SINGULARITY_EOS_HPP

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/eos/eos.hpp>

using singularity::EOS;
using singularity::TableSplit;

#if defined(__cplusplus)
extern "C" {
#endif

EOS *init_sg_IdealGas(const Real gamma_minus_1, const Real c_v);

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5
EOS *init_sg_SpinerDependsRhoT(const char *filename, const int matid,
                               const TableSplit split);
#endif

const Real get_sg_InternalEnergyFromDensityTemperature(EOS *eos, const Real density,
                                                       const Real temperature);

const Real get_sg_PressureFromDensityTemperature(EOS *eos, const Real density,
                                                 const Real temperature);

const Real get_sg_BulkModulusFromDensityTemperature(EOS *eos, const Real density,
                                                    const Real temperature);

const Real get_sg_SpecificHeatFromDensityTemperature(EOS *eos, const Real density,
                                                     const Real temperature);

const Real get_sg_MeanAtomicMass(EOS *eos);

const Real get_sg_MeanAtomicNumber(EOS *eos);

void modify_sg_ScaledEOS(EOS *eos, const Real scale);

void modify_sg_ShiftedEOS(EOS *eos, const Real shift);

void modify_sg_BilinearRampEOS(EOS *eos, const Real r0, const Real a, const Real b, const Real c);

void finalize_sg_eos(EOS *eos);

#if defined(__cplusplus)
}
#endif

#endif // SINGULARITY_EOS_HPP
