#include "fisher.hpp"

int main(int arc, char *argv[])
{
    Parameters parameters{};

	parameters.a = std::stod(argv[1]);
	parameters.b = std::stod(argv[2]);
	parameters.c = std::stod(argv[3]);
	parameters.biast = std::stod(argv[4]);
	parameters.mu_p = std::stod(argv[5]);
	parameters.mu_t = std::stod(argv[6]);
	parameters.max_mut_p = std::stod(argv[7]);
	parameters.max_mut_t = std::stod(argv[8]);
    parameters.init_p = std::stod(argv[9]);
    parameters.init_t = std::stod(argv[10]);
    parameters.only_positive = std::stoi(argv[11]);
    parameters.file_name = argv[12];

    Fisher f{parameters};

    f.run_simulation();
}
