#ifndef _actual_network_H_
#define _actual_network_H_

#include <network_properties.H>
#include <actual_network_data.H>

namespace Species {
    enum NetworkSpecies {He4 = 1, C12, O16, Ne20, Mg24, Si28, Ni56, NumberSpecies = Ni56};
};

const std::string network_name = "iso7";

void actual_network_init()
{
    actual_network_data_init();
}
