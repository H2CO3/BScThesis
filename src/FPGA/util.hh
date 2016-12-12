//
// util.hh
//
// General utility functions, types and macros
//
// Created on 07/03/2015
// by Arpad Goretity
//

#ifndef SWPARA_UTIL_HH
#define SWPARA_UTIL_HH

#include <cstddef>
#include <cstdlib>
#include <cstdint>

#ifndef __SYNTHESIS__
#include <cstdio>
#include <cassert>
#endif

#ifndef NDEBUG
#include <cstdarg>
#endif


#ifdef __SYNTHESIS__
#undef assert
#define assert(x)
#endif

#ifndef __SYNTHESIS__
#if __GNUC__
__attribute__((format(printf, 1, 2)))
#endif // __GNUC__
static inline void hls_debug(const char *fmt, ...)
{
#ifdef NDEBUG
	(void)fmt;
#else // NDEBUG
	std::va_list args;
	va_start(args, fmt);
	std::vfprintf(stderr, fmt, args);
	va_end(args);
#endif // NDEBUG
}
#else // __SYNTHESIS__
#define hls_debug(...)
#endif // __SYNTHESIS__

#define ARRAY_COUNT(a) (sizeof(a) / sizeof((a)[0]))

// TODO: is this good enough?
// XXX: this MUST be a signed integral type
typedef std::int16_t index_type;

// this MUST be an unsigned integer type,
// of which the size is the same as that of index_type.
typedef std::uint16_t size_type;

static_assert(sizeof(index_type) == sizeof(size_type), "size_type must be of the same size as index_type");

typedef std::uint32_t seq_count_type;


template<typename T, index_type N>
constexpr index_type array_count(const T (&)[N])
{
	return N;
}

#endif // SWPARA_UTIL_HH
