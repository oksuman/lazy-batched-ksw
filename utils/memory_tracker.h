#pragma once

#include <cstddef>
#include <sys/resource.h>

inline size_t getPeakRSSBytes() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    return ru.ru_maxrss;           // macOS: already in bytes
#else
    return ru.ru_maxrss * 1024;    // Linux: in KB
#endif
}

inline double getPeakRSSMB() {
    return getPeakRSSBytes() / (1024.0 * 1024.0);
}
