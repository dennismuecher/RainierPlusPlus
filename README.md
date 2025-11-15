# RAINIER++ - Modern C++ Edition

Complete rewrite of RAINIER in modern C++ with improved structure and maintainability. This is based on the original RAINIER.C by Leo Kirsch, which came in a single C source file. Leo Kirsch did a fantastic job and many people use the code to simulate level schemes and decays of nuclei via statistical models. However, the code was difficult to understand and maintain for other people, so here is an attempt to make RAINIER ready for the future.  

I used Claude Sonet 4.5 to translate RAINIER into C++. 

## Quick Start

1. **Prerequisites**
   - C++17 compiler
   - CMake 3.15+
   - ROOT 6.x

2. **Build**
   ```bash
   source /path/to/root/bin/thisroot.sh
   ./build.sh
   ```

3. **Run**
   ```bash
   cd build
   ./rainier
   ```

## Current Status

✅ **Phase 1 Complete**
- Modern C++ architecture
- Core classes implemented
- Build system working
- Configuration framework

🚧 **Phase 2 In Progress**
- Physics models (need implementation)
- Simulation engine (stub only)
- I/O handlers (basic stubs)

## Structure

```
rainier2/
├── include/          # Header files
│   ├── core/        # Level, Nucleus, Transition
│   ├── models/      # Physics models
│   ├── simulation/  # MC engine
│   └── io/          # Input/output
├── src/             # Implementation files
└── config/          # Configuration files
```

## What Works Now

- ✅ Compiles and runs
- ✅ Basic level loading (dummy data)
- ✅ Configuration system
- ✅ Clean class hierarchy

## What's Needed

- Physics model implementations
- Level file parsing
- Decay simulation
- ROOT output

See IMPLEMENTATION_GUIDE.md for details on completing each component.

## Documentation

- `README.md` - This file
- `PROJECT_STRUCTURE.md` - Complete file layout
- `IMPLEMENTATION_GUIDE.md` - How to implement remaining pieces
- `GETTING_STARTED.md` - Detailed startup guide

## Contact

For questions about implementation, see the documentation files or the inline code comments.
# RainierPlusPlus
