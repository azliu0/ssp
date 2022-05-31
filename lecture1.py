from temperature import *

cur_temp = get_cur_temp()


def F_to_C(temp_F):
    temp_C = (temp_F-32)*5/9
    return temp_C

def C_to_F(temp_C):
    temp_F = temp_C * 9/5 + 32
    return temp_F


cur_temps = [get_cur_temp1(),get_cur_temp2(),get_cur_temp3(),get_cur_temp4(),get_cur_temp5()]
print(cur_temps)

# for index in range(len(cur_temps)):
#     cur_temps[index] = C_to_F(cur_temps[index])

# seq_from_0_to_num = range(len(cur_temps))
# print(seq_from_0_to_num)
# print(list(seq_from_0_to_num))

cur_temps_F = []
for temp in cur_temps:
    cur_temps_F.append(C_to_F(temp)) # append is a method (method = function associated with a particular value)
print(cur_temps_F)

# temp = 87.8
# cur_temps_F.index(temp) #returns the first index where temp occurs in the list
# cur_temps_F.remove(temp) #removes the first instance of temp in the list (returns None)
# cur_temps_F.count(temp) #returns the number of occurrances of temp in the list
# cur_temps_F.sort() #modifies the list to be in ascending order (returns None)
# if temp in cur_temps_F:
#     print("found it!")


# print(cur_temps[-1])
# print(cur_temps[len(cur_temps)-1])



# hot_temp = 80
# hot_temp_C = (hot_temp - 32) * 5/9

# print("hot temp is", hot_temp_C, "degrees Celcius")

# a="hello"
# b="goodbye"
# print(a,b)

numbers = [0,1,2,3,4]
for i in range(1, len(numbers)):
    numbers[i] = numbers[i-1]
print(numbers)


while True:
    print(1)