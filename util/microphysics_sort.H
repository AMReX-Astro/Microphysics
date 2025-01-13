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
void swap(T& a, T& b) noexcept {
    // Helper function for sorting functions
    // Swaps the position of the two input

    T temp = a;
    a = b;
    b = temp;
}


template <typename T, int l, int m>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void bubbleSort_Array1D(amrex::Array1D<T, l, m>& sort_array,
                        const bool ascending=true) {

    // a simple bubble sort, for sorting Array1D
    // assumes index start from int l, and end at int m

    // int scratch;
    for (int i = l; i <= m; ++i) {
        for (int j = l; j <= m-i; ++j) {
            if (!compare(sort_array(j), sort_array(j+1), ascending)) {
                swap(sort_array(j), sort_array(j+1));
            }
        }
    }
}


template <typename T, typename P, int l, int m>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
int partition(amrex::Array1D<T, l, m>& sort_array,
              amrex::Array1D<P, l, m>& metric_array,
              int low, int high, const bool ascending=true) {
    // Helper function for quickSort
    // It uses the metric array as metric, but swaps sort_array

    // Choose pivot to be the last element

    P pivot = metric_array(high);
    int i = low - 1;
    for (int j = low; j <= high; j++) {
        if (compare(metric_array(j), pivot, ascending)) {
            i++;
            swap(sort_array(i), sort_array(j));
            swap(metric_array(i), metric_array(j));
        }
    }
    swap(sort_array(i+1), sort_array(high));
    swap(metric_array(i+1), metric_array(high));
    return i+1;
}


template <typename T, typename P, int l, int m>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void quickSort_Array1D(amrex::Array1D<T, l, m>& sort_array,
                       amrex::Array1D<P, l, m>& metric_array,
                       const bool ascending=true) {
    // quickSort implementation
    // This implementation uses the metric array as the sorting metric
    // to sort BOTH the sort_array and metric_array

    // Create a stack to keep track of the low and high index to the
    // left and right of the pivot

    amrex::Array1D<int, l, m> stack {0};
    int top = l - 1;

    // Set initial values of the range.

    stack(++top) = l;
    stack(++top) = m;

    // Keep popping from stack while it is not empty

    while (top >= l) {

        // Get high and low

        int high = stack(top--);
        int low = stack(top--);

        // Get the pivot index, which splits the array into two parts

        int pi = partition(sort_array, metric_array, low, high, ascending);

        // Record the low and high index to the left of the pivot

        if (pi - 1 > low) {
            stack(++top) = low;
            stack(++top) = pi - 1;
        }

        // Record the low and high index to the right of the pivot

        if (pi + 1 < high) {
            stack(++top) = pi + 1;
            stack(++top) = high;
        }
    }

}


#endif
