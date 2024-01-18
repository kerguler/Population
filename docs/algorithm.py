import numpy
from matplotlib import pyplot as plt

def briere1(T,T1,T2,a):
    v = a*T*(T-T1)*numpy.sqrt(T2-T)
    v[(T<T1) | (T>T2)] = 0.0
    return v

xr = numpy.linspace(0,40,1000)

plt.plot(xr,briere1(xr,10,30,0.01))
plt.show()

tm = numpy.array([20.0])
for i in range(10):
    tm = numpy.hstack([tm, tm[-1] + 0.2 - 0.01*tm[-1] + 1.0*numpy.sqrt((0.2+0.01*tm[-1]))*numpy.random.normal(scale=1.0)])

plt.plot(tm)
plt.show()

numpy.cumsum(briere1(tm,10,30,0.01))
numpy.round(numpy.cumsum(briere1(tm,10,30,0.01)))

plt.plot(numpy.cumsum(briere1(tm,10,30,0.01)))
plt.plot(numpy.round(numpy.cumsum(briere1(tm,10,30,0.01))))
plt.show()
