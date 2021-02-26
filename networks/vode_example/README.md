# vode_example

This is an implementation of the example stiff system from the VODE
paper / code comments.  Given an initial state Y = (y1, y2, y3) = (1,
0, 0), this has the long term behavior that y1, y2 = 0, y3 -> 1.

The VODE code comments give the following intermediate values:

```
   At t =  4.0000e-01   y =  9.851680e-01  3.386314e-05  1.479817e-02
   At t =  4.0000e+00   y =  9.055255e-01  2.240539e-05  9.445214e-02
   At t =  4.0000e+01   y =  7.158108e-01  9.184883e-06  2.841800e-01
   At t =  4.0000e+02   y =  4.505032e-01  3.222940e-06  5.494936e-01
   At t =  4.0000e+03   y =  1.832053e-01  8.942690e-07  8.167938e-01
   At t =  4.0000e+04   y =  3.898560e-02  1.621875e-07  9.610142e-01
   At t =  4.0000e+05   y =  4.935882e-03  1.984013e-08  9.950641e-01
   At t =  4.0000e+06   y =  5.166183e-04  2.067528e-09  9.994834e-01
   At t =  4.0000e+07   y =  5.201214e-05  2.080593e-10  9.999480e-01
   At t =  4.0000e+08   y =  5.213149e-06  2.085271e-11  9.999948e-01
   At t =  4.0000e+09   y =  5.183495e-07  2.073399e-12  9.999995e-01
   At t =  4.0000e+10   y =  5.450996e-08  2.180399e-13  9.999999e-01
```

