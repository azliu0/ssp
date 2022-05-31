#1.3 Exercise 1

def fib(n):
    if n==1: 
        return [0]
    if n==2:
        return [0,1]
    
    numbers = [0,1]
    fib1 = numbers[-2]
    fib2 = numbers[-1]
    for i in range(n-2):
        numbers.append(fib1+fib2)
        fib1 = numbers[-2]
        fib2 = numbers[-1]
    return numbers

def is_sorted1(numList):
    sortedList = numList
    sortedList.sort()
    if(sortedList == numList):
        return True
    return False







def is_sorted2(numList):
    if(numList ==  sorted(numList)):
        return True
    return False

print(is_sorted2([3,2,1]))
