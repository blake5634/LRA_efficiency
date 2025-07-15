import numpy as np

TRAPEZOIDAL = False

def integrate(x,xm1):
    if TRAPEZOIDAL:
        v =  0.5*(x+xm1)
        return v
    else:
        return x

if TRAPEZOIDAL:
    print('\nTrapezoidal Test\n')
else:
    print('\nRectangular Test\n')

print ('0,0 ', integrate(0,0))
print ('1,0 ', integrate(1,0))
print ('1,1 ' , integrate(1,1))
print ('1,-1', integrate(1,-1))


assert integrate(0,0)==0
if TRAPEZOIDAL:
    assert integrate(1,0)==0.5
    assert integrate(1,-1)==0
else:
    assert integrate(1,0)== 1
    assert integrate(1,-1)==1
assert integrate(1,1)==1

print('all tests passed')
