gen_point_pressure.py
=====================

- Calculates gravitational spherical harmonic coefficients for pressure values at individual points assuming a disc geometry

#### Calling Sequence
```python
from model_harmonics.gen_point_pressure import gen_point_pressure
Ylms = gen_point_pressure(P, G, R, lon, lat, AREA=AREA, LMAX=LMAX, LOVE=(hl,kl,ll))
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/model_harmonics/gen_point_pressure.py)

#### Inputs
 - `P`: Pressure (Pa)
 - `G`: Gravitational acceleration (m/s<sup>2</sup>)
 - `R`: Earth's radius at each data point (m)
 - `lon`: longitude of points
 - `lat`: latitude of points

#### Options
 - `AREA`: Area of each pressure cell (m<sup>2</sup>)
 - `LMAX`:  maximum spherical harmonic degree of the output harmonics
 - `MMAX`: maximum spherical harmonic order of the output harmonics
 - `LOVE`: input load Love numbers up to degree of truncation

#### Outputs
 - `clm`: Cosine spherical harmonic coefficients (geodesy normalization)
 - `slm`: Sine spherical harmonic coefficients (geodesy normalization)
 - `l`: spherical harmonic degree to LMAX
 - `m`: spherical harmonic order to MMAX
