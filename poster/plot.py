from pylab import *
import json as js

times = js.load(open('times.js'))
sizes = js.load(open('sizes.js'))

figure(1)
plot(
  sizes, times['cpp'], 
  sizes, times['for'], 
  sizes, times['pyt']
)
legend(['C++','FORTRAN','Python'])
show()
