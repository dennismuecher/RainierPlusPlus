# RAINIER 2.0 Project Structure

## Directory Layout

```
rainier2/
├── CMakeLists.txt        # Build configuration
├── build.sh              # Build script
├── README.md             # Main documentation
├── include/              # All headers
│   ├── Config.h
│   ├── core/            # Core classes
│   ├── models/          # Physics models
│   ├── simulation/      # Simulation engine
│   ├── io/              # Input/output
│   └── utils/           # Utilities
├── src/                  # Implementations
├── config/               # Configuration files
└── build/                # Build output (created)
```

## Implementation Status

### Complete
- Core level classes (Level, DiscreteLevel, ContinuumLevel)
- Transition class with selection rules
- Configuration framework
- Build system

### Stub/Partial
- Nucleus (basic structure, needs full implementation)
- Config (works but no JSON parsing yet)
- DecaySimulator (stub only)
- OutputManager (stub only)

### Not Yet Started
- LevelDensity models
- GammaStrength models
- SpinCutoff models
- InternalConversion
- Level file parsing
- Full simulation loop

## Next Steps

1. Implement LevelDensity.cpp
2. Implement level file reader
3. Implement simulation engine
4. Add physics models
5. Complete output system
