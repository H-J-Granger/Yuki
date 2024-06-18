#include <gmpxx.h>
#include <iostream>
#include <string>

int main() {
  std::ios_base::sync_with_stdio(false);
  for (int i = 0; i < 30; ++i) {
    std::string str, str2;
    std::cin >> str >> str2;
    mpz_class a(str, 16), b(str2, 16);
    mpz_class c = a * b;
    std::cout << c.get_str(16) << std::endl;
  }
  return 0;
}