gen_atmosphere_stokes.py
========================

 - Converts 3D atmospheric geopotential height and pressure difference fields from the spatial domain to spherical harmonic coefficients

#### Calling Sequence
```python
from model_harmonics.gen_atmosphere_stokes import gen_atmosphere_stokes
from gravity_toolkit.plm_holmes import plm_holmes
PLM,dPLM = plm_holmes(LMAX, np.cos(th))
Ylms = gen_atmosphere_stokes(GPH, pressure, lon, lat, LMAX=LMAX,
    ELLIPSOID='WGS84', GEOID=geoid, PLM=PLM, LOVE=(hl,kl,ll))
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/model_harmonics/gen_atmosphere_stokes.py)

#### Arguments
- `GPH`: geopotential heights at model levels
- `pressure`: pressure differences between model levels
- `lon`: longitude array
- `lat`: latitude array

#### Keyword arguments
- `LMAX`:  maximum spherical harmonic degree of the output harmonics
- `MMAX`: maximum spherical harmonic order of the output harmonics
- `ELLIPSOID`: reference ellipsoid name
- `GEOID`: geoid height
- `PLM`: input Legendre polynomials (for improving computational time)
- `LOVE`: input load Love numbers up to degree of truncation
- `METHOD`: method of integrating over pressure levels
    * `'SW02'`: [Swenson and Wahr (2002)](https://doi.org/10.1029/2000JB000024)
    * `'BC05'`: [Boy and Chao (2005)](https://doi.org/10.1029/2002JB002333)

#### Returns
- `clm`: Cosine spherical harmonic coefficients (geodesy normalization)
- `slm`: Sine spherical harmonic coefficients (geodesy normalization)
- `l`: spherical harmonic degree to LMAX
- `m`: spherical harmonic order to MMAX

