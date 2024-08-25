#include <actual_network.H>


namespace network
{

    AMREX_GPU_MANAGED Array1D<Real, 1, 10> semenov_x;
    AMREX_GPU_MANAGED Array1D<Real, 1, 1000> semenov_y;
    AMREX_GPU_MANAGED Array2D<Real, 1, 10, 1, 1000> semenov_z;

}

void actual_network_init()
{

	using namespace network;
    std::string filename = "Semenov_PlanckOpacity.dat";
    std::cout << "Reading tables from " << filename << std::endl;

    // Open the file and check if it exists
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("ERROR: file " + filename + " not found!");
    }

    std::string line;
    
    // Skip comments and read until the first non-comment line
    while (std::getline(file, line)) {
        if (line[0] != '#') {
            break;
        }
    }

    // Check if the first line is formatted correctly
    if (line.find(',') == std::string::npos) {
        throw std::runtime_error("ERROR: file " + filename + " should contain the number of rows and columns in the format 'RR, CC'");
    }

    // Process data line by line
    for (int i = 1; i <= semenov_x.size(); ++i) {
        for (int j = 1; j <= semenov_y.size(); ++j) {
            if (!std::getline(file, line)) {
                throw std::runtime_error("ERROR: Unexpected end of file while reading data.");
            }

            // Trim the line to remove any leading/trailing whitespace
            line.erase(0, line.find_first_not_of(" \t"));  // Remove leading whitespace
            line.erase(line.find_last_not_of(" \t") + 1);  // Remove trailing whitespace

            std::istringstream iss(line);
            Real x_val, y_val, z_val;
            if (!(iss >> x_val >> y_val >> z_val)) {
                std::cout << "line: " << line << std::endl;
                throw std::runtime_error("ERROR: Insufficient data on line " + std::to_string(i * semenov_y.size() + j + 1));
            }

            semenov_x(i) = x_val;
            semenov_y(j) = y_val;
            semenov_z(i,j) = z_val;
        }

    }


}

