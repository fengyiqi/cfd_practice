#include <fstream>
#include "output.h"

OutputWriter::OutputWriter(const Block& block, std::string file_name) : block_(block), file_name_(file_name) {}


void OutputWriter::WriteCSV() {
    std::ofstream csv_file;
    csv_file.open(file_name_, std::ios::out);
    csv_file << "x,density,velocityx,pressure,energy" << std::endl;
    for (unsigned int i = 0; i < GI::TCX(); i++){
        csv_file << block_.x_coordinate_buffer_[i] << ",";
        csv_file << block_.primitive_buffer_[PIndex(FI::PrimeStateEnum::Density)][i] << ",";
        csv_file << block_.primitive_buffer_[PIndex(FI::PrimeStateEnum::VelocityX)][i] << ",";
        csv_file << block_.primitive_buffer_[PIndex(FI::PrimeStateEnum::Pressure)][i] << ",";
        csv_file << block_.conservative_buffer_[EIndex(FI::EquationEnum::Energy)][i] << std::endl;
    }
    csv_file.close();
}