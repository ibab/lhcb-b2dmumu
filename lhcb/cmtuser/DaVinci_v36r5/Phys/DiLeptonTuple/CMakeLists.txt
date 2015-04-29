################################################################################
# Package: DiLeptonTuple
################################################################################
gaudi_subdir(DiLeptonTuple v1r4)

gaudi_depends_on_subdirs(Phys/DaVinciMCKernel
                         Phys/DecayTreeTupleBase)

gaudi_add_module(DiLeptonTuple
                 src/*.cpp
                 LINK_LIBRARIES DaVinciMCKernelLib DecayTreeTupleBaseLib)

gaudi_install_python_modules()

