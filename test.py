import random, time;

start = time.time()
for i in range(30000):
  a = int(input(), 16)
  b = int(input(), 16)
  c = a * b
  print(hex(c))
end = time.time()