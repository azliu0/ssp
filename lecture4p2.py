vec = [5,5,5]
q = vec # makes q another name (alias) for the list
q[1] = 7
print(vec) #so vec changes


foo = [4,4,4]
bar = foo
print("id(foo) =", id(foo))
print("id(bar) =", id(bar))
foo = [1,2,3]

print(bar) # only changes foo, does not change bar
print("id(foo) =", id(foo))
print("id(bar) =", id(bar))

a=5
b=5
print("id(a) =", id(a))
print("id(b) =", id(b))
b = b+1
print("id(a) =", id(a))
print("id(b) =", id(b))
b=b-1
a=a+1
print("id(a) =", id(a))
print("id(b) =", id(b))
