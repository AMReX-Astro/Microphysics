#ifndef UNIT_TEST_H
#define UNIT_TEST_H

void write_job_info(const std::string& dir);

extern "C" {

void runtime_pretty_print(int* jobinfo_file_name, const int* jobinfo_file_length);

}

#endif
