from . import molwt

M_ATMOS   = 5.1352e18 # mass of atmosphere, kg

# Conversion between ppm CO2 and GtC emissions
ppm_gtc   = M_ATMOS/1e18*molwt.C/molwt.AIR
