#ifndef YUKI_LIB_FIX_TYPES_H_
#define YUKI_LIB_FIX_TYPES_H_

#include <cstdint>
#include <cstdlib>

namespace yuki {

using ::std::size_t;
using uint64 = ::std::uint64_t;
using int64 = ::std::int64_t;

#ifndef _GLIBCXX_USE_INT128
#error Yuki library needs 128-bit integer type.
#else
using uint128 = __uint128_t;
#endif

template <typename T>
class type_display;

}  // namespace yuki

#endif  // #ifndef YUKI_LIB_FIX_TYPES_H_