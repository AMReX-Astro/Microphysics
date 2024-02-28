#ifndef MICROPHYSICS_SORT_H
#define MICROPHYSICS_SORT_H

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool compare(const T& a, const T& b, const bool ascending=true) {
    // Helper function for sorting functions
    // compares two values
    // asecending for least to greatest

    if (ascending) {
        return a < b;
    }
    else {
        return a > b;
    }
}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void swap(T& a, T& b) {
    // Helper function for sorting functions
    // Swaps the position of the two input

    T temp = a;
    a = b;
    b = temp;
}

template <typename T, int l, int m>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void bubbleSort_Array1D(amrex::Array1D<T, l, m>& array_1d,
                        const bool ascending=true) {

    // a simple bubble sort, for sorting Array1D
    // assumes index start from int l, and end at int m

    // int scratch;
    for (int i = l; i <= m; ++i) {
        for (int j = l; j <= m-i; ++j) {
            if (!compare(array_1d(j), array_1d(j+1), ascending)) {
                swap(array_1d(j), array_1d(j+1));
                // scratch = array_1d(j);
                // array_1d(j) = array_1d(j+1);
                // array_1d(j+1) = scratch;
            }
        }
    }
}

template <typename T, int l, int m>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
int partition(amrex::Array1D<T, l, m>& array_1d,
              const int low, const int high, const bool ascending=true) {
    // Helper function for quickSort

    T pivot = array_1d(high);
    int i = low - 1;
    for (int j = low; j <= high; j++) {
        if (compare(array_1d(j), pivot, ascending)) {
            i++;
            swap(array_1d(i), array_1d(j));
        }
    }
    swap(array_1d(i+1), array_1d(high));
    return i+1;
}

template <typename T, int l, int m>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void quickSort_Array1D(amrex::Array1D<T, l, m>& array_1d,
                       const int low, const int high,
                       const bool ascending=true) {
    // quickSort implementation
    // This version uses the array itself as the metric for sorting
    if (low < high) {

        // pi is the partition return index of pivot

        int pi = partition(array_1d, low, high);

        // Recursive quick sort calls.
        quickSort_Array1D(array_1d, low, pi-1, ascending);
        quickSort_Array1D(array_1d, pi+1, high, ascending);
    }
}


// template <int l, int m>
// AMREX_GPU_HOST_DEVICE AMREX_INLINE
// void quickSort_Array1D(amrex::Array1D<int, l, m>& array_1d,
//                        amrex::Array1D<int, l, m>& metric_array) {
//     // This version uses the metric_array as the metric for sorting
// }



#endif
