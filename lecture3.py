import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

def C_to_F(temp_C):
    temp_F = temp_C * 9/5 + 32
    return temp_F


nums = np.array([6,8,12,22])
nums_slice = nums[1:3]
print(nums_slice)
print(nums)
# print(nums.toList())
print(nums)
print(nums*3)
print(nums + [1,1,0,0])
print(nums - 3)
print(3-nums)
print(C_to_F(nums))
print(np.zeros(100))


mat = np.array([[ 4, 5, 6],[ 7, 8, 9],[15,16,17]] )

print(mat)
print(mat*3)
print(mat+3)


mat2 = np.arange(10,22).reshape(3,4)
print(mat2)
print(mat2.shape)
print(mat2.shape[0])
print(mat2.shape[1])
print(mat2.size)
print(mat2.T)
print(mat2 @ mat2.T)

print(np.array([1,2,3]))




""" restrictions on arrays
while lists can grow (via append, +, etc), arrays are fixed in size
arrays can only contain one type of data (all numbers, or all strings, or all booleans)

"""


# #LIST[start, end, step]
# #default: [0, past the end, 1], when the step is positive
# #default: [end, 0, -1], when the step is negative

# nums = [6,8,12,22]

# print(nums * 3) #creates 3 copies in one list
# print(nums + [1,1,0,0]) #appends [1,1,0,0] to the end of the list


# nums[:] #entire list
# nums[1]

# nums[1:3:1] #[6,8,12]
# nums[1:3] #[6,8,12] default step = 1
# nums[:3] #[6,8,12], default start = 1


# nums[-4:-2] #[6,8]

# nums[::-1] # [22,12,8,6] reversed the whole list


# mat = [[ 4, 5, 6], # 0
#        [ 7, 8, 9], # 1
#        [15,16,17]] # 2

# mat[1][0] # index ROW, then COLUMN

# print(mat[1:])




happy_file = open("happy_little_file.txt", "r")
# happy_lines = happy_file.readlines()
# print(happy_lines)

i=0
for line in happy_file:
    print("line", i, "is", line.strip())
    i = i+1


# happy_text = happy_file.read()
# print(happy_text)
# "\n" is the character for a new line


np.savetxt("data.txt", nums)
nums_from_file=np.loadtxt("data.txt")
print(nums_from_file)

new_file = open("another_happy_file.txt", "w")
new_file.write("just a happy little file\n")
new_file.write("life as a file is good\n")
new_file.close()

# from astropy.io import fits
# import matplotlib.pyplot as plt

image = fits.getdata("SSP2021Astro_Test_Image.FIT")
print(image)

print(image.shape)
print(image.size)
print(image.max())
print(image.mean())

# plt.imshow(image)
# plt.gray()
# plt.show()

# box = image[285:295, 305:315]
# print(box)
# plt.imshow(box)
# plt.show()