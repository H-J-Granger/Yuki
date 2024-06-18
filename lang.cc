#include <algorithm>
#include <cmath>
#include <cstdio>

template <size_t len1, size_t len2>
requires(len1 < len2) void func(void) {
  std::puts("a");
}
template <size_t len1, size_t len2>
requires(len1 >= len2) void func(void) {
  std::puts("b");
}

int main() {
  func<2, 3>();
  func<3, 2>();
  return 0;
}