#include <chrono>
#include <iostream>
#include <random>
#include <string>

std::mt19937 gen(std::chrono::duration_cast<std::chrono::nanoseconds>(
                     std::chrono::steady_clock::now().time_since_epoch())
                     .count());

std::string produce_hex(int len, bool sign = true) {
  std::string ret;
  ret.resize(len + sign);
  int n = gen() % len + 1;
  for (int i = 0; i < n; ++i) {
    int x = gen() % 16;
    ret[i] = (x < 10 ? '0' : 'a' - 10) + x;
  }
  return (sign && gen() % 2 ? "-" : "") + ret.substr(0, n);
}

int main() {
  std::ios_base::sync_with_stdio(false);
  for (int i = 0; i < 30; ++i) {
    std::cout << produce_hex(50000) << std::endl;
    std::cout << produce_hex(50000) << std::endl;
  }
  return 0;
}