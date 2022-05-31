from decimal import Decimal



x=0.1
print(Decimal("0.1"))
print(Decimal("100000000000000000000000000.1"))
print(format(x,".30f"))

print(0.1+0.1+0.1 == 0.3)
print(format(0.1+0.1+0.1, ".30f"))
print(format(0.3, ".30f"))
print(format(1000.1, ".30f"))
print(format(100000000000.1, ".30f"))
print(format(100000000000000000000000000.1, ".30f"))

a = (2.5*0.1)*1.5
b = (2.5)*(0.1*1.5)
print(a==b)
threshold = 1e-9
print(abs(a-b) < threshold)
print(format(a, ".30f"))
print(format(b, ".30f"))
