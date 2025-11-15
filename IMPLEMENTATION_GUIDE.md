# RAINIER 2.0 Implementation Guide

## Quick Reference

### Files That Need Implementation

1. **Priority 1 - Core Functionality**
   - `src/core/Nucleus.cpp` - Level loading (partially done)
   - `src/io/LevelFileReader.cpp` - Parse level files (NEW)
   - `src/models/LevelDensity.cpp` - BSFG/CTM models (NEW)

2. **Priority 2 - Physics**
   - `src/models/SpinCutoff.cpp` - Spin distributions (NEW)
   - `src/models/GammaStrength.cpp` - GSF models (NEW)

3. **Priority 3 - Simulation**
   - `src/simulation/DecaySimulator.cpp` - Full implementation
   - `src/simulation/CascadeEvent.cpp` - Event storage (NEW)

4. **Priority 4 - I/O**
   - `src/Config.cpp` - Add JSON parsing
   - `src/io/OutputManager.cpp` - ROOT histograms

## How to Add a New Component

1. Create header in `include/appropriate_dir/`
2. Create implementation in `src/appropriate_dir/`
3. Add both to CMakeLists.txt
4. Rebuild: `cd build && cmake .. && make`

## Testing Individual Components

You can test components independently:

```cpp
// test_level.cpp
#include "core/DiscreteLevel.h"
int main() {
    DiscreteLevel level(1.5, 2.0, 1, 1e6);
    std::cout << "Energy: " << level.getEnergy() << "\n";
    return 0;
}
```

Compile:
```bash
g++ -std=c++17 -I include test_level.cpp src/core/Level.cpp \
    -o test_level `root-config --cflags --libs`
```

## Reference Original RAINIER

For physics formulas, see the original RAINIER.C:
- Lines 159-230: Level density models
- Lines 232-295: Spin cutoff models
- Lines 363-472: Gamma strength functions
- Lines 503-600: Width calculations
- Lines 1212-1450: Cascade simulation

## Getting Help

Each header file has detailed documentation. The implementations should follow the patterns shown in `src/core/Level.cpp`.
