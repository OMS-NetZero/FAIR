from . import molwt

M_ATMOS          = 5.1352e18 # mass of atmosphere, kg
EARTH_RADIUS     = 6371000   # m
SECONDS_PER_YEAR = 60 * 60 * 24 * 365.24219 # Length of tropical year
    # https://en.wikipedia.org/wiki/Tropical_year

# Conversion between ppm CO2 and GtC emissions
ppm_gtc   = M_ATMOS/1e18*molwt.C/molwt.AIR
