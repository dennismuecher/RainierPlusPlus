// TransitionPool.h - Memory pool for reusing Transition objects
#ifndef RAINIER_TRANSITION_POOL_H
#define RAINIER_TRANSITION_POOL_H

#include "core/Transition.h"
#include "core/ContinuumLevel.h"
#include <vector>
#include <memory>

namespace rainier {

/**
 * @brief Object pool for Transition objects to avoid allocations in hot loops
 * 
 * Pre-allocates a pool of Transition objects and reuses them instead of
 * constantly allocating/deallocating. This eliminates the massive allocation
 * overhead seen in profiling (33% of runtime in operator new).
 * 
 * Usage:
 *   TransitionPool pool(1000);  // Pre-allocate 1000 transitions
 *   
 *   // In hot loop:
 *   Transition* trans = pool.acquire();  // Reuse existing object
 *   trans->setInitialLevel(level1);
 *   trans->setFinalLevel(level2);
 *   // ... use transition ...
 *   
 *   pool.reset();  // Ready for next iteration
 */
class TransitionPool {
public:
    /**
     * @brief Construct pool with initial capacity
     * @param initialSize Number of transitions to pre-allocate
     */
    explicit TransitionPool(size_t initialSize = 1000) {
        pool_.reserve(initialSize);
        
        // Create dummy levels for pool objects (will be overwritten when used)
        // These are stored as members to avoid recreating them during growth
        dummyInit_ = std::make_shared<ContinuumLevel>(0.0, 0.0, 1, 0, 0);
        dummyFinal_ = std::make_shared<ContinuumLevel>(0.0, 0.0, 1, 0, 0);
        
        for (size_t i = 0; i < initialSize; ++i) {
            // Use new + unique_ptr constructor to avoid make_unique template issues
            pool_.push_back(std::unique_ptr<Transition>(
                new Transition(dummyInit_, dummyFinal_, 0.0, 0.0)
            ));
        }
        nextIndex_ = 0;
    }
    
    /**
     * @brief Get a transition from the pool (reuses existing if available)
     */
    Transition* acquire() {
        if (nextIndex_ >= pool_.size()) {
            // Need to grow pool - allocate in chunks for efficiency
            size_t oldSize = pool_.size();
            size_t newSize = oldSize + std::max(size_t(100), oldSize / 2);
            pool_.reserve(newSize);
            
            for (size_t i = oldSize; i < newSize; ++i) {
                pool_.push_back(std::unique_ptr<Transition>(
                    new Transition(dummyInit_, dummyFinal_, 0.0, 0.0)
                ));
            }
        }
        
        return pool_[nextIndex_++].get();
    }
    
    /**
     * @brief Reset pool for reuse (call after processing all transitions)
     */
    void reset() {
        nextIndex_ = 0;
    }
    
    /**
     * @brief Get current pool size
     */
    size_t size() const { return pool_.size(); }
    
    /**
     * @brief Get number of transitions currently in use
     */
    size_t inUse() const { return nextIndex_; }

private:
    std::vector<std::unique_ptr<Transition>> pool_;
    size_t nextIndex_;
    
    // Dummy levels used for pool initialization
    // These are never actually used - real levels set via setInitialLevel/setFinalLevel
    std::shared_ptr<ContinuumLevel> dummyInit_;
    std::shared_ptr<ContinuumLevel> dummyFinal_;
};

} // namespace rainier

#endif // RAINIER_TRANSITION_POOL_H
