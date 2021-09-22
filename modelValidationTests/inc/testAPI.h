//
// Created by mbarbone on 7/22/21.
//

#ifndef DPM_MODELVALIDATIONAPI_TESTAPI_H_
#define DPM_MODELVALIDATIONAPI_TESTAPI_H_

#include <string>
#include <tuple>
#include <utility>
#include <vector>

class AbstractTest {
   protected:
    const std::string gInputDataDir;
    const double gPrimaryEnergy;
    const int gNumPrimaries;
    const int gMaterialIndex;
    // histogram parameter
    int hbins;
    double xmin;
    double xmax;
    double hdel;
    double ihdel;
    std::vector<double> theHist;
    std::vector<double> theCosineHist;

   public:
    explicit AbstractTest(std::string gInputDataDir = "../data",
                          double gPrimaryEnergy     = 12.345,
                          int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0,
                          int hbins = 101);
    void writeHists(const std::string &filenamePrefix);
    __attribute__((unused)) std::vector<std::tuple<double, double>>
    getEnergyHist();
    __attribute__((unused)) std::vector<std::tuple<double, double>>
    getCosHist();
    virtual void simulate() = 0;
};

class BremTest : public AbstractTest {
   public:
    explicit BremTest(std::string gInputDataDir = "../data",
                      double gPrimaryEnergy     = 12.345,
                      int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0,
                      int hbins = 101)
        : AbstractTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries,
                       gMaterialIndex, hbins) {}
    void simulate() final;
    void writeHists();
};

class MollerTest : public AbstractTest {
   public:
    explicit MollerTest(std::string gInputDataDir = "../data",
                        double gPrimaryEnergy     = 12.345,
                        int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0,
                        int hbins = 101)
        : AbstractTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries,
                       gMaterialIndex, hbins) {}
    void simulate() final;
    void writeHists();
};

class MSCAngularDeflectionTest : public AbstractTest {
   public:
    explicit MSCAngularDeflectionTest(std::string gInputDataDir = "../data",
                                      double gPrimaryEnergy     = 12.345,
                                      int gNumPrimaries         = 1.0E+5,
                                      int gMaterialIndex = 0, int hbins = 101)
        : AbstractTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries,
                       gMaterialIndex, hbins) {}
    void simulate() final;
    void writeHists();
};

class AnnihilationTest : public AbstractTest {
   public:
    explicit AnnihilationTest(std::string gInputDataDir = "../data",
                              double gPrimaryEnergy     = 12.345,
                              int gNumPrimaries         = 1.0E+5,
                              int gMaterialIndex = 0, int hbins = 101)
        : AbstractTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries,
                       gMaterialIndex, hbins) {}
    void simulate() final;
    void writeHists();
};

class ComptonPrimaryTest : public AbstractTest {
   public:
    explicit ComptonPrimaryTest(std::string gInputDataDir = "../data",
                                double gPrimaryEnergy     = 12.345,
                                int gNumPrimaries         = 1.0E+5,
                                int gMaterialIndex = 0, int hbins = 101)
        : AbstractTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries,
                       gMaterialIndex, hbins) {}
    void simulate() final;
    void writeHists();
    ~ComptonPrimaryTest() = default;
};

class ComptonSecondaryTest : public AbstractTest {
   public:
    explicit ComptonSecondaryTest(std::string gInputDataDir = "../data",
                                  double gPrimaryEnergy     = 12.345,
                                  int gNumPrimaries         = 1.0E+5,
                                  int gMaterialIndex = 0, int hbins = 101)
        : AbstractTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries,
                       gMaterialIndex, hbins) {}
    void simulate() final;
    void writeHists();
    ~ComptonSecondaryTest() = default;
};

class SimElectronData;
class SimMaterialData;
class Geom;

class ElectronTrackingTest {
   private:
    SimMaterialData *theMaterialData;
    SimElectronData *elData;
    Geom *geom;

   public:
    ElectronTrackingTest(const std::string &g_input_data_dir, double voxelSize);
    virtual ~ElectronTrackingTest();

   public:
    __attribute__((unused)) int Traking(const double *position,
                                        const double *direction, double eKin,
                                        double &numTr1MFP, double &numMollerMFP,
                                        double invMollerMFP,
                                        double &numBremMFP);
};

#endif  // DPM_MODELVALIDATIONAPI_TESTAPI_H_
