An implementation of the analytical EOS for hypervelocity impacts
provided in Tillotson (1962):

https://ui.adsabs.harvard.edu/abs/1962geat.rept.3216T/abstract

This EOS primarily serves to give P = P(rho, e) for giant impact
simulations. Since the Tillotson formulation can sometimes yield
negative pressure, we make sure to apply a floor in that case.

The implementation here follows Reinhardt and Stadel (2017):

https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.4252R/abstract

Particularly, several quantities (particularly sound speed and
temperature) are provided either in Appendix A of that paper and
this GitHub repo linked in the paper:

https://github.com/chreinhardt/tillotson

The default parameters for this EOS are currently set to the values
for granite from Reinhardt and Stadel's Table A1. At present there is
no support for multi-material compositions.
