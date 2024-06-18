#include "lib/fix_int.h"

#include <iostream>
#include <string>
#include <vector>

yuki::fix_int<32 * 6400> a, b;
int main() {
  std::ios_base::sync_with_stdio(false);
  std::string str;
  for (int i = 0; i < 30; ++i) {
    std::cin >> str;
    a.from_hex_string(str);
    std::cin >> str;
    b.from_hex_string(str);
    // int k;
    // std::cin >> k;
    auto c = a * b;
    std::cout << c.to_hex_string() << std::endl;
  }
  std::cout << std::endl;
  return 0;
}
