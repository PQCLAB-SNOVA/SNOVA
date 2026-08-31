/**
 * @file ct_poison.h
 */
#ifndef SNOVA_CT_POISON_H
#define SNOVA_CT_POISON_H

#if defined(SNOVA_CT_TEST)

#if defined(__has_feature)
#  if __has_feature(memory_sanitizer)
#    define SNOVA_CT_MSAN 1
#  endif
#endif

#if defined(SNOVA_CT_REQUIRE_MSAN) && !defined(SNOVA_CT_MSAN)
#error "SNOVA_CT_REQUIRE_MSAN set but MemorySanitizer is NOT enabled (missing -fsanitize=memory?) -- refusing to build a vacuous CT gate"
#endif

#if defined(SNOVA_CT_MSAN)
#include <sanitizer/msan_interface.h>
#define SNOVA_CT_POISON(p, n) __msan_poison((p), (n))
#define SNOVA_CT_DECLASSIFY(p, n) __msan_unpoison((p), (n))
#define SNOVA_CT_ASSERT_PUBLIC_MEM(p, n) __msan_check_mem_is_initialized((p), (n))
#else
#include <valgrind/memcheck.h>
#define SNOVA_CT_POISON(p, n) VALGRIND_MAKE_MEM_UNDEFINED((p), (n))
#define SNOVA_CT_DECLASSIFY(p, n) VALGRIND_MAKE_MEM_DEFINED((p), (n))
#define SNOVA_CT_ASSERT_PUBLIC_MEM(p, n) ((void)VALGRIND_CHECK_MEM_IS_DEFINED((p), (n)))
#endif

#else
#define SNOVA_CT_POISON(p, n) do { } while (0)
#define SNOVA_CT_DECLASSIFY(p, n) do { } while (0)
#define SNOVA_CT_ASSERT_PUBLIC_MEM(p, n) do { } while (0)
#endif

#endif
