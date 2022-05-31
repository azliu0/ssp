import math

def getChange(change):
    numQuarters = change // 25
    change = change % 25
    numDimes = change // 10 
    change = change % 10 
    numNickels = change // 5 
    change = change % 5
    numPennies = change
    print("Quarters:", numQuarters, "Dimes:", numDimes, "Nickels:", numNickels, "Pennies:", numPennies)

getChange(141)

def distance(x1, y1, x2,y2):
    return math.sqrt((x1-x2)**2+(y1-y2)**2)
print(distance(0,0,3,4))
    