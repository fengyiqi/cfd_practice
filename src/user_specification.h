#ifndef USER_SPECIFICATION_H
#define USER_SPECIFICATION_H

#include <array>

class GridInformation {

    static constexpr unsigned int internal_cells = 3000;
    static constexpr unsigned int ghost_cells = 4;
    static constexpr unsigned int total_cells = 2 * ghost_cells + internal_cells;

public:
    GridInformation() = delete;
    ~GridInformation() = default;
    GridInformation(GridInformation const& ) = delete;
    GridInformation& operator=(GridInformation const&) = delete;
    GridInformation(GridInformation &&) = delete;
    GridInformation& operator=(GridInformation&&) = delete;
    // return internal cells number
    static constexpr unsigned int ICX() { return internal_cells; }
    // return ghost cells number
    static constexpr unsigned int GCX() { return ghost_cells; }
    // return total cells (including internal cells and two-side ghost cells)
    static constexpr unsigned int TCX() { return total_cells; }
    // return the index of the first internal cell
    static constexpr unsigned int FICX() { return ghost_cells; }
    // return the index of the last internal cell
    static constexpr unsigned int LIXC() { return ghost_cells + internal_cells - 1; }
    // return the index of the first ghost cell on the righ side
    static constexpr unsigned int FRGX() { return ghost_cells + internal_cells; }

};

using GI = GridInformation;

namespace FieldInformation{
    enum class States {
        Primitives,
        Conservatives
    };
    // define the valid equations
    enum class EquationEum : unsigned int{
        Mass,
        MomentumX,
        Energy,
    };
    // define the valid primetive states
    enum class PrimeStateEum : unsigned int{
        Density,
        VelocityX,
        Pressure,
    };
    // assemble all valid equations into an array for range loop
    constexpr std::array<std::underlying_type<EquationEum>::type, 3> Equations{
        static_cast<std::underlying_type<EquationEum>::type>(EquationEum::Mass),
        static_cast<std::underlying_type<EquationEum>::type>(EquationEum::MomentumX),
        static_cast<std::underlying_type<EquationEum>::type>(EquationEum::Energy)
    };
    // assemble all valid prime states into an array for range loop
    constexpr std::array<std::underlying_type<PrimeStateEum>::type, 3> PrimeStates{
        static_cast<std::underlying_type<PrimeStateEum>::type>(PrimeStateEum::Density),
        static_cast<std::underlying_type<PrimeStateEum>::type>(PrimeStateEum::VelocityX),
        static_cast<std::underlying_type<PrimeStateEum>::type>(PrimeStateEum::Pressure)
    };
    // get the numbers
    constexpr unsigned int EN() { return Equations.size(); }
    constexpr unsigned int PN() { return PrimeStates.size(); }

}

namespace FI = FieldInformation;

constexpr std::underlying_type<FI::EquationEum>::type EIndex(const FI::EquationEum e) {
    return static_cast<std::underlying_type<FI::EquationEum>::type>( e );
}
constexpr std::underlying_type<FI::PrimeStateEum>::type PIndex(const FI::PrimeStateEum p) {
    return static_cast<std::underlying_type<FI::PrimeStateEum>::type>( p );
}

#endif