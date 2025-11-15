# Getting Started with RAINIER 2.0

## What You Have

A complete C++ framework with:
- ✅ Working build system
- ✅ Core classes implemented
- ✅ Clean architecture
- ✅ Ready for physics implementation

## Build and Run

```bash
# Source ROOT
source /path/to/root/bin/thisroot.sh

# Build
./build.sh

# Run
cd build
./rainier
```

## Expected Output

You should see:
```
╔════════════════════════════════════════════════════════╗
║   RAINIER 2.0 - Modern C++ Edition                    ║
╚════════════════════════════════════════════════════════╝

Loading configuration...
=== Initializing Nucleus ===
Creating dummy discrete levels...
=== Starting Simulations ===
[stub implementations run]
=== Simulation Complete ===
```

## What's Next?

The framework is ready. Now we need to implement:
1. Physics models (level density, gamma strength)
2. Level file parsing
3. Simulation engine
4. Output system

See IMPLEMENTATION_GUIDE.md for details.

## Customization

Edit `config/example.json` to change:
- Nucleus (Z, A)
- Number of realizations
- Events per realization
- Output files

## Help

- Check header files for API documentation
- See original RAINIER.C for physics
- Review example implementations in src/core/Level.cpp
