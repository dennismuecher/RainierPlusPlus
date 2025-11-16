#ifndef RAINIER_TRANSITION_H
#define RAINIER_TRANSITION_H

#include <memory>

namespace rainier {

class Level;

class Transition {
public:
    enum class Type {
        NONE = 0,
        E1 = 1,
        M1_E2 = 2,
        M1_PURE = 3,
        E2_PURE = 4
    };

    Transition(std::shared_ptr<Level> initialLevel,
               std::shared_ptr<Level> finalLevel,
               double branchingRatio = 1.0,
               double icc = 0.0);

    std::shared_ptr<Level> getInitialLevel() const { return initialLevel_; }
    std::shared_ptr<Level> getFinalLevel() const { return finalLevel_; }
    double getGammaEnergy() const;
    double getBranchingRatio() const { return branchingRatio_; }
    double getInternalConversionCoeff() const { return icc_; }
    Type getType() const { return type_; }
    double getMixingRatio() const { return mixingRatio_; }
    double getPartialWidth() const { return partialWidth_; }

    void setBranchingRatio(double br) { branchingRatio_ = br; }
    void setInternalConversionCoeff(double icc) { icc_ = icc; }
    void setType(Type type) { type_ = type; }
    void setMixingRatio(double delta) { mixingRatio_ = delta; }
    void setPartialWidth(double width) { partialWidth_ = width; }
	
	void setInitialLevel(std::shared_ptr<Level> level) { initialLevel_ = level; }
	void setFinalLevel(std::shared_ptr<Level> level) { finalLevel_ = level; }



    static Type determineType(double spinI, int parityI,
                             double spinF, int parityF,
                             bool isEvenA);

    bool isAllowed() const { return type_ != Type::NONE; }
    double conversionProbability() const;

private:
    std::shared_ptr<Level> initialLevel_;
    std::shared_ptr<Level> finalLevel_;
    double branchingRatio_;
    double icc_;
    Type type_;
    double mixingRatio_;
    double partialWidth_;
};

} // namespace rainier

#endif
