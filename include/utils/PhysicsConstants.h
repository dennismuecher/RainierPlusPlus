#ifndef RAINIER_PHYSICS_CONSTANTS_H
#define RAINIER_PHYSICS_CONSTANTS_H

namespace rainier {
namespace constants {

constexpr double PI = 3.14159265358979323846;
constexpr double FOUR_PI_SQUARED = 4.0 * PI * PI;
constexpr double E_EULER = 2.71828182845904523536;

constexpr double HBAR = 6.5821195e-7;  // MeV * fs
constexpr double HBAR_C = 197.327;      // MeV * fm

constexpr double K_X1 = 1.0 / (3.0 * PI * PI * HBAR * HBAR * 197.327 * 197.327);
constexpr double K_X2 = 5.204155555e-08;

constexpr double FS_TO_PS = 1e-3;
constexpr double FS_TO_NS = 1e-6;
constexpr double LOG_2 = 0.693147180559945309417;

constexpr double DEFAULT_MAX_HALFLIFE = 1e9;
constexpr double MIN_ENERGY = 1e-8;
constexpr double MIN_DENSITY = 1e-10;

} // namespace constants
} // namespace rainier

#endif
