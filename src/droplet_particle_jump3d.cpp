/*
 * OpenLB case: two droplets coalescing/jumping on a superhydrophobic wall while
 * carrying one spherical hydrophilic solid particle.
 *
 * The fluid part follows the OpenLB free-energy two-component examples. The
 * particle is a compact Lagrangian closure that samples the OpenLB phase and
 * velocity fields and writes a VTK-compatible trajectory/diagnostic stream.
 */

#include "olb3D.h"
#include "olb3D.hh"

#include <array>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace olb;
using namespace olb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19<>

namespace {

struct Vec3 {
  T x = 0.;
  T y = 0.;
  T z = 0.;

  Vec3 operator+(const Vec3& rhs) const { return {x + rhs.x, y + rhs.y, z + rhs.z}; }
  Vec3 operator-(const Vec3& rhs) const { return {x - rhs.x, y - rhs.y, z - rhs.z}; }
  Vec3 operator*(T scale) const { return {x * scale, y * scale, z * scale}; }
  Vec3 operator/(T scale) const { return {x / scale, y / scale, z / scale}; }
};

T dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
T norm(const Vec3& v) { return std::sqrt(dot(v, v)); }
Vec3 normalize(const Vec3& v)
{
  const T n = norm(v);
  return n > 1.e-12 ? v / n : Vec3{0., 0., 0.};
}

struct CaseConfig {
  int nx = 160;
  int ny = 96;
  int nz = 96;
  int resolution = 160;
  bool periodicX = true;
  bool periodicY = true;

  T dropletRadius = 22.;
  Vec3 drop1 = {58., 48., 21.};
  Vec3 drop2 = {102., 48., 21.};
  T interfaceWidth = 4.;
  T initialBridge = 1.5;
  T liquidDensity = 1.;
  T gasDensity = 0.02;

  T substrateContactAngleDeg = 165.;
  int wallMaterial = 2;

  T particleRadius = 8.;
  Vec3 particleCenter = {80., 48., 23.};
  T particleDensityRatio = 1.8;
  T particleContactAngleDeg = 35.;
  T dragCoefficient = 0.35;
  T capillaryAdhesion = 0.18;
  T restitution = 0.15;

  T alpha = 1.;
  T kappa1 = 0.005;
  T kappa2 = 0.005;
  T gamma = 10.;
  T relaxationTime = 1.;

  int maxIter = 70000;
  int vtkIter = 500;
  int statIter = 100;
  int csvIter = 20;
  std::string outputDir = "./tmp/droplet_particle_jump3d/";
};

std::string trim(std::string value)
{
  const auto begin = value.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) {
    return "";
  }
  const auto end = value.find_last_not_of(" \t\r\n");
  return value.substr(begin, end - begin + 1);
}

bool parseBool(const std::string& value)
{
  return value == "true" || value == "1" || value == "yes" || value == "on";
}

CaseConfig readConfig(const std::string& path)
{
  CaseConfig c;
  if (path.empty()) {
    return c;
  }

  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("Could not open config file: " + path);
  }

  std::string section;
  std::string line;
  while (std::getline(input, line)) {
    const auto comment = line.find('#');
    if (comment != std::string::npos) {
      line = line.substr(0, comment);
    }
    line = trim(line);
    if (line.empty()) {
      continue;
    }
    if (line.front() == '[' && line.back() == ']') {
      section = trim(line.substr(1, line.size() - 2));
      continue;
    }

    const auto equals = line.find('=');
    if (equals == std::string::npos) {
      continue;
    }
    const std::string key = trim(line.substr(0, equals));
    const std::string value = trim(line.substr(equals + 1));
    const std::string full = section + "." + key;
    const T scalar = std::atof(value.c_str());
    const int integer = std::atoi(value.c_str());

    if (full == "domain.nx") c.nx = integer;
    else if (full == "domain.ny") c.ny = integer;
    else if (full == "domain.nz") c.nz = integer;
    else if (full == "domain.resolution") c.resolution = integer;
    else if (full == "domain.periodicX") c.periodicX = parseBool(value);
    else if (full == "domain.periodicY") c.periodicY = parseBool(value);
    else if (full == "droplets.radius") c.dropletRadius = scalar;
    else if (full == "droplets.centerX1") c.drop1.x = scalar;
    else if (full == "droplets.centerX2") c.drop2.x = scalar;
    else if (full == "droplets.centerY") c.drop1.y = c.drop2.y = scalar;
    else if (full == "droplets.centerZ") c.drop1.z = c.drop2.z = scalar;
    else if (full == "droplets.interfaceWidth") c.interfaceWidth = scalar;
    else if (full == "droplets.initialBridge") c.initialBridge = scalar;
    else if (full == "droplets.liquidDensity") c.liquidDensity = scalar;
    else if (full == "droplets.gasDensity") c.gasDensity = scalar;
    else if (full == "substrate.contactAngleDeg") c.substrateContactAngleDeg = scalar;
    else if (full == "substrate.wallMaterial") c.wallMaterial = integer;
    else if (full == "particle.radius") c.particleRadius = scalar;
    else if (full == "particle.centerX") c.particleCenter.x = scalar;
    else if (full == "particle.centerY") c.particleCenter.y = scalar;
    else if (full == "particle.centerZ") c.particleCenter.z = scalar;
    else if (full == "particle.densityRatio") c.particleDensityRatio = scalar;
    else if (full == "particle.contactAngleDeg") c.particleContactAngleDeg = scalar;
    else if (full == "particle.dragCoefficient") c.dragCoefficient = scalar;
    else if (full == "particle.capillaryAdhesion") c.capillaryAdhesion = scalar;
    else if (full == "particle.restitution") c.restitution = scalar;
    else if (full == "freeEnergy.alpha") c.alpha = scalar;
    else if (full == "freeEnergy.kappa1") c.kappa1 = scalar;
    else if (full == "freeEnergy.kappa2") c.kappa2 = scalar;
    else if (full == "freeEnergy.gamma") c.gamma = scalar;
    else if (full == "freeEnergy.relaxationTime") c.relaxationTime = scalar;
    else if (full == "run.maxIter") c.maxIter = integer;
    else if (full == "run.vtkIter") c.vtkIter = integer;
    else if (full == "run.statIter") c.statIter = integer;
    else if (full == "run.csvIter") c.csvIter = integer;
    else if (full == "run.outputDir") c.outputDir = value;
  }
  return c;
}

T smoothSphere(const Vec3& p, const Vec3& center, T radius, T width)
{
  const T signedDistance = radius - norm(p - center);
  return 0.5 * (1. + std::tanh(signedDistance / std::max(width, T(1.e-6))));
}

T liquidIndicator(const CaseConfig& c, const Vec3& p)
{
  const T d1 = smoothSphere(p, c.drop1, c.dropletRadius, c.interfaceWidth);
  const T d2 = smoothSphere(p, c.drop2, c.dropletRadius, c.interfaceWidth);
  const Vec3 bridgeCenter = (c.drop1 + c.drop2) * 0.5;
  const T bridge = smoothSphere(p, bridgeCenter, c.initialBridge, c.interfaceWidth);
  return std::min(T(1.), d1 + d2 + bridge);
}

class InitialRho final : public AnalyticalF3D<T,T> {
public:
  explicit InitialRho(const CaseConfig& config) : AnalyticalF3D<T,T>(1), _c(config) { }
  bool operator()(T output[], const T input[]) override
  {
    const T chi = liquidIndicator(_c, {input[0], input[1], input[2]});
    output[0] = _c.gasDensity + (_c.liquidDensity - _c.gasDensity) * chi;
    return true;
  }
private:
  const CaseConfig& _c;
};

class InitialPhi final : public AnalyticalF3D<T,T> {
public:
  explicit InitialPhi(const CaseConfig& config) : AnalyticalF3D<T,T>(1), _c(config) { }
  bool operator()(T output[], const T input[]) override
  {
    const T chi = liquidIndicator(_c, {input[0], input[1], input[2]});
    output[0] = 2. * chi - 1.;
    return true;
  }
private:
  const CaseConfig& _c;
};

struct ParticleState {
  Vec3 position;
  Vec3 velocity;
  T radius = 1.;
  T mass = 1.;
};

class ParticleCoupler {
public:
  explicit ParticleCoupler(const CaseConfig& c)
    : _c(c)
  {
    _particle.position = c.particleCenter;
    _particle.radius = c.particleRadius;
    _particle.mass = (4./3.) * M_PI * std::pow(c.particleRadius, 3) * c.particleDensityRatio;
  }

  const ParticleState& particle() const { return _particle; }

  void advance(const Vec3& liquidCom, const Vec3& liquidVelocity, T localPhi)
  {
    const T liquidContact = std::max(T(0.), std::min(T(1.), 0.5 * (localPhi + 1.)));
    const Vec3 drag = (liquidVelocity - _particle.velocity) * (_c.dragCoefficient * liquidContact);
    const Vec3 adhesionDir = normalize(liquidCom - _particle.position);
    const T hydrophilicGain = std::max(T(0.), std::cos(_c.particleContactAngleDeg * M_PI / 180.));
    const Vec3 capillary = adhesionDir * (_c.capillaryAdhesion * hydrophilicGain * liquidContact / _particle.mass);

    _particle.velocity = _particle.velocity + drag + capillary;
    _particle.position = _particle.position + _particle.velocity;

    const T wallLimit = _particle.radius + 0.5;
    if (_particle.position.z < wallLimit) {
      _particle.position.z = wallLimit;
      if (_particle.velocity.z < 0.) {
        _particle.velocity.z = -_c.restitution * _particle.velocity.z;
      }
    }
  }

private:
  const CaseConfig& _c;
  ParticleState _particle;
};

void prepareGeometry(SuperGeometry3D<T>& superGeometry, const CaseConfig& c, UnitConverter<T,DESCRIPTOR>& converter)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare geometry" << std::endl;

  superGeometry.rename(0, 2);
  Vector<T,3> extend(c.nx + 2., c.ny + 2., c.nz - converter.getPhysDeltaX());
  Vector<T,3> origin(-1., -1., 0.5 * converter.getPhysDeltaX());
  IndicatorCuboid3D<T> fluidDomain(extend, origin);
  superGeometry.rename(2, 1, fluidDomain);
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();
}

T wettingParameterFromContactAngle(T contactAngleDeg, T alpha, T kappa)
{
  const T centered = std::max(T(-0.98), std::min(T(0.98), std::cos(contactAngleDeg * M_PI / 180.)));
  return -0.25 * alpha * kappa * centered;
}

void prepareLattice(
  SuperLattice3D<T,DESCRIPTOR>& sLattice1,
  SuperLattice3D<T,DESCRIPTOR>& sLattice2,
  Dynamics<T,DESCRIPTOR>& bulkDynamics1,
  Dynamics<T,DESCRIPTOR>& bulkDynamics2,
  UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>& superGeometry,
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sOnBC1,
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sOnBC2,
  const CaseConfig& c)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare lattice" << std::endl;

  sLattice1.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T,DESCRIPTOR>());
  sLattice2.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T,DESCRIPTOR>());
  sLattice1.defineDynamics(superGeometry, 1, &bulkDynamics1);
  sLattice2.defineDynamics(superGeometry, 1, &bulkDynamics2);
  sLattice1.defineDynamics(superGeometry, c.wallMaterial, &instances::getNoDynamics<T,DESCRIPTOR>());
  sLattice2.defineDynamics(superGeometry, c.wallMaterial, &instances::getNoDynamics<T,DESCRIPTOR>());

  const T h1 = wettingParameterFromContactAngle(c.substrateContactAngleDeg, c.alpha, c.kappa1);
  const T h2 = -h1;
  sOnBC1.addFreeEnergyWallBoundary(superGeometry, c.wallMaterial, c.alpha, c.kappa1, c.kappa2, h1, h2, 1);
  sOnBC2.addFreeEnergyWallBoundary(superGeometry, c.wallMaterial, c.alpha, c.kappa1, c.kappa2, h1, h2, 2);

  std::vector<T> zero(3, T());
  AnalyticalConst3D<T,T> zeroVelocity(zero);
  InitialRho rho(c);
  InitialPhi phi(c);

  sLattice1.iniEquilibrium(superGeometry, 1, rho, zeroVelocity);
  sLattice2.iniEquilibrium(superGeometry, 1, phi, zeroVelocity);
  sLattice1.iniEquilibrium(superGeometry, c.wallMaterial, rho, zeroVelocity);
  sLattice2.iniEquilibrium(superGeometry, c.wallMaterial, phi, zeroVelocity);
  sLattice1.initialize();
  sLattice2.initialize();
  sLattice1.communicate();
  sLattice2.communicate();
}

void prepareCoupling(SuperLattice3D<T,DESCRIPTOR>& sLattice1, SuperLattice3D<T,DESCRIPTOR>& sLattice2,
                     SuperGeometry3D<T>& superGeometry, const CaseConfig& c)
{
  FreeEnergyChemicalPotentialGenerator3D<T,DESCRIPTOR> chemical(c.alpha, c.kappa1, c.kappa2);
  FreeEnergyForceGenerator3D<T,DESCRIPTOR> force;
  chemical.shift(0, 0, 0);
  force.shift(0, 0, 0);
  sLattice1.addLatticeCoupling(superGeometry, 1, chemical, sLattice2);
  sLattice2.addLatticeCoupling(superGeometry, 1, force, sLattice1);
}

struct LiquidMoments {
  Vec3 com;
  T volume = 0.;
  bool touchesWall = true;
};

LiquidMoments sampleLiquidMoments(SuperLattice3D<T,DESCRIPTOR>& sLattice2, const CaseConfig& c)
{
  SuperLatticeDensity3D<T,DESCRIPTOR> phiField(sLattice2);
  AnalyticalFfromSuperF3D<T,T> phi(phiField, true, 1);
  Vec3 weighted{0., 0., 0.};
  T volume = 0.;
  bool wallContact = false;
  for (int ix = 0; ix < c.nx; ix += 2) {
    for (int iy = 0; iy < c.ny; iy += 2) {
      for (int iz = 0; iz < c.nz; iz += 2) {
        T input[3] = {T(ix), T(iy), T(iz)};
        T value[1] = {0.};
        phi(value, input);
        const T chi = std::max(T(0.), std::min(T(1.), 0.5 * (value[0] + 1.)));
        const Vec3 p{T(ix), T(iy), T(iz)};
        weighted = weighted + p * chi;
        volume += chi;
        wallContact = wallContact || (iz <= 2 && chi > 0.25);
      }
    }
  }
  LiquidMoments moments;
  moments.com = volume > 0. ? weighted / volume : Vec3{T(c.nx)/2., T(c.ny)/2., T(c.nz)/2.};
  moments.volume = volume * 8.;
  moments.touchesWall = wallContact;
  return moments;
}

Vec3 sampleFluidVelocity(SuperLattice3D<T,DESCRIPTOR>& sLattice1, const Vec3& position)
{
  SuperLatticeVelocity3D<T,DESCRIPTOR> velocityField(sLattice1);
  AnalyticalFfromSuperF3D<T,T> velocity(velocityField, true, 1);
  T input[3] = {position.x, position.y, position.z};
  T output[3] = {0., 0., 0.};
  velocity(output, input);
  return {output[0], output[1], output[2]};
}

T samplePhi(SuperLattice3D<T,DESCRIPTOR>& sLattice2, const Vec3& position)
{
  SuperLatticeDensity3D<T,DESCRIPTOR> phiField(sLattice2);
  AnalyticalFfromSuperF3D<T,T> phi(phiField, true, 1);
  T input[3] = {position.x, position.y, position.z};
  T output[1] = {0.};
  phi(output, input);
  return output[0];
}

void writeParticleVtk(const CaseConfig& c, const ParticleState& particle, int iT)
{
  std::ostringstream name;
  name << c.outputDir << "particle_" << iT << ".vtp";
  std::ofstream out(name.str());
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "<PolyData><Piece NumberOfPoints=\"1\" NumberOfVerts=\"1\">\n";
  out << "<PointData Scalars=\"radius\"><DataArray type=\"Float64\" Name=\"radius\" format=\"ascii\">"
      << particle.radius << "</DataArray></PointData>\n";
  out << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
      << particle.position.x << " " << particle.position.y << " " << particle.position.z
      << "</DataArray></Points>\n";
  out << "<Verts><DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">0</DataArray>"
      << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">1</DataArray></Verts>\n";
  out << "</Piece></PolyData></VTKFile>\n";
}

class DiagnosticsWriter {
public:
  explicit DiagnosticsWriter(const std::string& outputDir)
    : _csv(outputDir + "diagnostics.csv")
  {
    _csv << "step,time,liquid_com_x,liquid_com_y,liquid_com_z,liquid_volume,"
         << "particle_x,particle_y,particle_z,particle_vx,particle_vy,particle_vz,takeoff\n";
  }

  void write(int iT, T time, const Vec3& liquidCom, T liquidVolume, const ParticleState& particle, bool takeoff)
  {
    _csv << iT << ',' << time << ','
         << liquidCom.x << ',' << liquidCom.y << ',' << liquidCom.z << ',' << liquidVolume << ','
         << particle.position.x << ',' << particle.position.y << ',' << particle.position.z << ','
         << particle.velocity.x << ',' << particle.velocity.y << ',' << particle.velocity.z << ','
         << (takeoff ? 1 : 0) << '\n';
  }
private:
  std::ofstream _csv;
};

void getResults(
  SuperLattice3D<T,DESCRIPTOR>& sLattice1,
  SuperLattice3D<T,DESCRIPTOR>& sLattice2,
  int iT,
  SuperGeometry3D<T>& superGeometry,
  Timer<T>& timer,
  UnitConverter<T,DESCRIPTOR>& converter,
  ParticleCoupler& particleCoupler,
  DiagnosticsWriter& diagnostics,
  const CaseConfig& c)
{
  SuperVTMwriter3D<T> vtmWriter("droplet_particle_jump3d");
  if (iT == 0) {
    SuperLatticeGeometry3D<T,DESCRIPTOR> geometry(sLattice1, superGeometry);
    SuperLatticeCuboid3D<T,DESCRIPTOR> cuboid(sLattice1);
    SuperLatticeRank3D<T,DESCRIPTOR> rank(sLattice1);
    vtmWriter.write(geometry);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }

  if (iT % c.statIter == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice1.getStatistics().print(iT, converter.getPhysTime(iT));
    sLattice2.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  if (iT % c.vtkIter == 0) {
    SuperLatticeVelocity3D<T,DESCRIPTOR> velocity(sLattice1);
    SuperLatticeDensity3D<T,DESCRIPTOR> rho(sLattice1);
    SuperLatticeDensity3D<T,DESCRIPTOR> phi(sLattice2);
    rho.getName() = "rho";
    phi.getName() = "phi";
    velocity.getName() = "velocity";
    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(rho);
    vtmWriter.addFunctor(phi);
    vtmWriter.write(iT);
    writeParticleVtk(c, particleCoupler.particle(), iT);
  }

  if (iT % c.csvIter == 0) {
    const LiquidMoments moments = sampleLiquidMoments(sLattice2, c);
    const bool takeoff = !moments.touchesWall && moments.com.z > c.dropletRadius + 2.;
    diagnostics.write(iT, converter.getPhysTime(iT), moments.com, moments.volume, particleCoupler.particle(), takeoff);
  }
}

} // namespace

int main(int argc, char* argv[])
{
  olbInit(&argc, &argv);
  const std::string configPath = argc > 1 ? argv[1] : "configs/default.ini";
  const CaseConfig c = readConfig(configPath);
  singleton::directories().setOutputDir(c.outputDir);

  OstreamManager clout(std::cout, "main");
  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    (T)c.resolution,
    c.relaxationTime,
    (T)c.nx,
    0.0001,
    1.002e-8,
    1.0);
  converter.print();

  Vector<T,3> extend((T)c.nx, (T)c.ny, (T)c.nz);
  Vector<T,3> origin(0., 0., 0.);
  IndicatorCuboid3D<T> cuboid(extend, origin);
#ifdef PARALLEL_MODE_MPI
  CuboidGeometry3D<T> cGeometry(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());
#else
  CuboidGeometry3D<T> cGeometry(cuboid, converter.getPhysDeltaX());
#endif
  cGeometry.setPeriodicity(c.periodicX, c.periodicY, false);

  HeuristicLoadBalancer<T> loadBalancer(cGeometry);
  loadBalancer.print();
  SuperGeometry3D<T> superGeometry(cGeometry, loadBalancer);
  prepareGeometry(superGeometry, c, converter);

  SuperLattice3D<T,DESCRIPTOR> sLattice1(superGeometry);
  SuperLattice3D<T,DESCRIPTOR> sLattice2(superGeometry);
  ForcedBGKdynamics<T,DESCRIPTOR> bulkDynamics1(converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T,DESCRIPTOR>());
  FreeEnergyBGKdynamics<T,DESCRIPTOR> bulkDynamics2(converter.getLatticeRelaxationFrequency(), c.gamma, instances::getBulkMomenta<T,DESCRIPTOR>());

  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sOnBC1(sLattice1);
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sOnBC2(sLattice2);
  createLocalBoundaryCondition3D<T,DESCRIPTOR>(sOnBC1);
  createLocalBoundaryCondition3D<T,DESCRIPTOR>(sOnBC2);
  prepareLattice(sLattice1, sLattice2, bulkDynamics1, bulkDynamics2, converter, superGeometry, sOnBC1, sOnBC2, c);
  prepareCoupling(sLattice1, sLattice2, superGeometry, c);

  SuperExternal3D<T,DESCRIPTOR> sExternal1(superGeometry, sLattice1, sLattice1.getOverlap());
  SuperExternal3D<T,DESCRIPTOR> sExternal2(superGeometry, sLattice2, sLattice2.getOverlap());
  ParticleCoupler particleCoupler(c);
  DiagnosticsWriter diagnostics(c.outputDir);

  Timer<T> timer(c.maxIter, superGeometry.getStatistics().getNvoxel());
  timer.start();
  for (int iT = 0; iT <= c.maxIter; ++iT) {
    const LiquidMoments moments = sampleLiquidMoments(sLattice2, c);
    const Vec3 localVelocity = sampleFluidVelocity(sLattice1, particleCoupler.particle().position);
    const T localPhi = samplePhi(sLattice2, particleCoupler.particle().position);
    particleCoupler.advance(moments.com, localVelocity, localPhi);

    getResults(sLattice1, sLattice2, iT, superGeometry, timer, converter, particleCoupler, diagnostics, c);

    sLattice1.collideAndStream();
    sLattice2.collideAndStream();
    sLattice1.communicate();
    sLattice2.communicate();
    sLattice1.executeCoupling();
    sExternal1.communicate();
    sExternal2.communicate();
    sLattice2.executeCoupling();
  }

  timer.stop();
  timer.printSummary();
  return 0;
}
