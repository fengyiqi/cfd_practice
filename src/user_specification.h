#ifndef USER_SPECIFICATION_H
#define USER_SPECIFICATION_H

class GridInformation {

    static constexpr unsigned int internal_cells = 500;
    static constexpr unsigned int ghost_cells = 4;
    static constexpr unsigned int total_cells = 2 * ghost_cells + internal_cells;

public:
    GridInformation() = delete;
    ~GridInformation() = default;
    GridInformation(GridInformation const& ) = delete;
    GridInformation& operator=(GridInformation const&) = delete;
    GridInformation(GridInformation &&) = delete;
    GridInformation& operator=(GridInformation&&) = delete;

    static constexpr unsigned int ICX() { return internal_cells; }
    static constexpr unsigned int GCX() { return ghost_cells; }

    static constexpr unsigned int TCX() { return total_cells; }

    static constexpr unsigned int FICX() { return ghost_cells; }

    static constexpr unsigned int LIXC() { return ghost_cells + internal_cells - 1; }

    static constexpr unsigned int FHHX() { return ghost_cells + internal_cells; }

};

using GI = GridInformation;

enum States {
    Primitives,
    Conservatives
};

enum ConservativePool {
    Mass,
    MomentumX,
    Energy,
};

enum PrimeStatePool {
   Density,
   VelocityX,
   Pressure,
};

class FieldInformation {
    // Primitive states number
    static constexpr unsigned int primitive_states_number = 3;
    static constexpr unsigned int conservative_states_number = 3;

public:
    FieldInformation() = delete;
    ~FieldInformation() = default;
    FieldInformation(FieldInformation const& ) = delete;
    FieldInformation& operator=(FieldInformation const&) = delete;
    FieldInformation(FieldInformation &&) = delete;
    FieldInformation& operator=(FieldInformation&&) = delete;

    static constexpr unsigned int PN() { return primitive_states_number; }
    static constexpr unsigned int CN() { return conservative_states_number; }

};

using FI = FieldInformation;

#endif