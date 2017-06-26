# Import required packages
import numpy as np
from scipy.optimize import root

# Define a function which gives the relationship between iIRF_100 and scaling factor, alpha
def iirf100_interp_funct(alpha,a,tau,targ_iirf100):
    iirf100_arr = alpha*(np.sum(a*tau*(1.0 - np.exp(-100.0/(tau*alpha)))))
    return iirf100_arr   -  targ_iirf100

# Define the FAIR simple climate model function
def fair_scm(tstep=1.0,
             emissions=False,
             other_rf=0.0,
             co2_concs=False,
             q=np.array([0.33,0.41]),
             tcrecs=np.array([1.6,2.75]),
             d=np.array([239.0,4.1]),
             a=np.array([0.2173,0.2240,0.2824,0.2763]),
             tau=np.array([1000000,394.4,36.54,4.304]),
             r0=32.40,
             rC=0.019,
             rT=4.165,
             F_2x=3.74,
             C_0=278.0,
             ppm_gtc=2.123,
             iirf100_max=97.0,
             in_state=[[0.0,0.0,0.0,0.0],[0.0,0.0],0.0],
             restart_out=False):

  # # # ------------ CALCULATE Q ARRAY ------------ # # #
  # If TCR and ECS are supplied, overwrite the q array
  k = 1.0 - (d/70.0)*(1.0 - np.exp(-70.0/d))
  if type(tcrecs) in [np.ndarray,list]:
    q =  (1.0 / F_2x) * (1.0/(k[0]-k[1])) \
        * np.array([tcrecs[0]-k[1]*tcrecs[1],k[0]*tcrecs[1]-tcrecs[0]])

  # # # ------------ SET UP OUTPUT TIMESERIES VARIABLES ------------ # # #
  # the integ_len variable is used to store the length of our timeseries
  # by default FAIR is not concentration driven
  conc_driven=False
  # here we check if FAIR is emissions driven
  if type(emissions) in [np.ndarray,list]:
    integ_len = len(emissions)
    if (type(other_rf) in [np.ndarray,list]) and (len(other_rf)!=integ_len):
        raise ValueError("The emissions and other_rf timeseries don't have the same length")
    elif type(other_rf) in [int,float]:
        other_rf = np.full(integ_len,other_rf)
  
  # here we check if FAIR is concentration driven
  elif type(co2_concs) in [np.ndarray,list]:
    integ_len = len(co2_concs)
    conc_driven = True
    if (type(other_rf) in [np.ndarray,list]) and (len(other_rf)!=integ_len):
        raise ValueError("The concentrations and other_rf timeseries don't have the same length")
    elif type(other_rf) in [int,float]:
        other_rf = np.full(integ_len,other_rf)

  # finally we check if only a non-CO2 radiative forcing timeseries has been supplied
  elif type(other_rf) in [np.ndarray,list]:
    integ_len = len(other_rf)
    if type(emissions) in [int,float]:
        emissions = np.full(integ_len,emissions)
    else:
        emissions = np.zeros(integ_len)

  else:
    raise ValueError("Neither emissions, co2_concs or other_rf is defined as a timeseries")

  RF = np.zeros(integ_len)
  C_acc = np.zeros(integ_len)
  iirf100 = np.zeros(integ_len)

  carbon_boxes_shape = (integ_len,4)
  R_i = np.zeros(carbon_boxes_shape)
  C = np.zeros(integ_len)

  thermal_boxes_shape = (integ_len,2)
  T_j = np.zeros(thermal_boxes_shape)
  T = np.zeros(integ_len)

  # # # ------------ FIRST TIMESTEP ------------ # # #
  R_i_pre = in_state[0]
  C_pre = np.sum(R_i_pre) + C_0
  T_j_pre = in_state[1]
  C_acc_pre = in_state[2]

  if conc_driven:
    C[0] = co2_concs[0]
  
  else:
    # Calculate the parametrised iIRF and check if it is over the maximum 
    # allowed value
    iirf100[0] = r0 + rC*C_acc_pre + rT*np.sum(T_j_pre)
    if iirf100[0] >= iirf100_max:
      iirf100[0] = iirf100_max
      
    # Determine a solution for alpha using scipy's root finder
    time_scale_sf = (root(iirf100_interp_funct,0.16,args=(a,tau,iirf100[0])))['x']

    # Multiply default timescales by scale factor
    tau_new = time_scale_sf * tau

    # Compute the updated concentrations box anomalies from the decay of the 
    # previous year and the emisisons
    R_i[0] = R_i_pre*np.exp(-tstep/tau_new) \
              + (emissions[0,np.newaxis])*a*tau_new*(1-np.exp(-tstep/tau_new)) / ppm_gtc

    C[0] = np.sum(R_i[0]) + C_0

    # Calculate the additional carbon uptake
    C_acc[0] =  C_acc_pre + emissions[0] - (C[0]-(np.sum(R_i_pre) + C_0)) * ppm_gtc

  # Calculate the radiative forcing using the previous timestep's CO2 concentration

  RF[0] = (F_2x/np.log(2.)) * np.log(C_pre/C_0) + other_rf[0]

  # Update the thermal response boxes
  T_j[0] = RF[0,np.newaxis]*q*(1-np.exp((-tstep)/d)) + T_j_pre*np.exp(-tstep/d)

  # Sum the thermal response boxes to get the total temperature anomlay
  T[0] = np.sum(T_j[0])

  # # # ------------ REST OF RUN ------------ # # #
  for x in range(1,integ_len):
    if conc_driven:
      C[x] = co2_concs[x]
    
    else:
      # Calculate the parametrised iIRF and check if it is over the maximum 
      # allowed value
      iirf100[x] = r0 + rC*C_acc[x-1] + rT*T[x-1]
      if iirf100[x] >= iirf100_max:
        iirf100[x] = iirf100_max
        
      # Determine a solution for alpha using scipy's root finder
      time_scale_sf = (root(iirf100_interp_funct,time_scale_sf,args=(a,tau,iirf100[x])))['x']

      # Multiply default timescales by scale factor
      tau_new = time_scale_sf * tau

      # Compute the updated concentrations box anomalies from the decay of the previous year and the emisisons
      R_i[x] = R_i[x-1]*np.exp(-tstep/tau_new) \
              + (emissions[x,np.newaxis])*a*tau_new*(1-np.exp(-tstep/tau_new)) / ppm_gtc

      # Sum the boxes to get the total concentration anomaly
      C[x] = np.sum(R_i[x]) + C_0

      # Calculate the additional carbon uptake
      C_acc[x] =  C_acc[x-1] + emissions[x] * tstep - (C[x]-C[x-1]) * ppm_gtc

    # Calculate the radiative forcing using the previous timestep's CO2 concentration
    RF[x] = (F_2x/np.log(2.)) * np.log((C[x-1]) /C_0) + other_rf[x]

    # Update the thermal response boxes
    T_j[x] = T_j[x-1]*np.exp(-tstep/d) + RF[x,np.newaxis]*q*(1-np.exp(-tstep/d))
    
    # Sum the thermal response boxes to get the total temperature anomaly
    T[x] = np.sum(T_j[x])

  if restart_out:
    return C, T, (R_i[-1],T_j[-1],C_acc[-1])
  else:
    return C, T

# Define a function to plot FAIR's inputs and outputs
# Shifts variables to appropriately represent their definitions (continuous 
# throughout timestep or end of timestep)
def plot_fair(emms,
              conc,
              forc,
              temp,
              y_0=0,
              tuts=False,
              infig=False,
              ax1in=None,
              ax2in=None,
              ax3in=None,
              ax4in=None,
              colour={'emms':'black',
                     'conc':'blue',
                     'forc':'orange',
                     'temp':'red'},
              label=None,
              linestyle='-',
             ):
  """
  One line summary 

  More details of behaviour if required.

  # # ------------ ARGUMENTS ------------ # #
  # sublime snippet for variable description is 'vardesc'

  * => Optional argument
  ^ => Keyword argument

  # # ------------ RETURN VALUE ------------ # #
  # sublime snippet for variable description is 'vardesc'
  # fig: (matplotlib axes object)
  #   the matplotlib axes with all our plotted variables

  # # ------------ SIDE EFFECTS ------------ # #
  # document side effects here

  # # ------------ EXCEPTIONS ------------ # #
  # sublime snippet for exception description is 'excdesc'

  # # ------------ RESTRICTIONS ------------ # #
  Document any restrictions on when the function can be called

  """

  # One line break before anything else

  # # # ------------ IMPORT REQUIRED MODULES ------------ # # #
  # # ------------ STANDARD LIBRARY ------------ # #
  from math import ceil

  # # ------------ THIRD PARTY ------------ # #
  import numpy as np

  from matplotlib import pyplot as plt
  plt.style.use('seaborn-darkgrid')
  plt.rcParams['figure.figsize'] = 16, 9
  plt.rcParams['lines.linewidth'] = 1.5

  font = {'weight' : 'normal',
          'size'   : 12}

  plt.rc('font', **font)

  # # ------------ LOCAL APPLICATION/LIBRARY SPECIFIC ------------ # #

  # # # ------------ CODE ------------ # # #
  # # ------------ SORT OUT INPUT VARIABLES ------------ # #
  pts = {'emms':emms,
         'forc':forc,
         'conc':conc,
         'temp':temp}
  
  integ_len = 0

  for j,var in enumerate(pts):
    if type(pts[var]) == list:
      pts[var] = np.array(pts[var])
      integ_len = len(pts[var])
    elif type(pts[var]) == np.ndarray:
      integ_len = len(pts[var])

  if integ_len == 0:
    for name in pts:
      print "{0}: {1}".format(name,type(pts[name]))
    raise ValueError("Error: I can't work out which one of your input variables is a timeseries")

  for j,var in enumerate(pts):
    if type(pts[var]) == np.ndarray and len(pts[var]) == integ_len:
      pass
    elif type(pts[var]) == np.ndarray and len(pts[var]) != integ_len:
      for name in pts:
        print "{0}: {1}\nlength: {2}".format(name,
                                             type(pts[name]),
                                             len(pts[name]))
      raise ValueError("Error: Your timeseries are not all the same length, I don't know what to do")
    else:
      if type(pts[var]) in [float,int]:
        pts[var] = np.full(integ_len,pts[var])
      else:
        pts[var] = np.zeros(integ_len)

  # # ------------ SORT OUT TIME VARIABLE ------------ # #
  # state variables are valid at the end of the timestep so we
  # go from 1 - integ_len + 1 rather than 0 - integ_len
  time = np.arange(0.99,integ_len+0.99) + y_0

  if not tuts:
    tuts = 'units unknown'

  # # ------------ PREPARE FLUX VARIABLES FOR PLOTTING  ------------ # #
  # Flux variables are assumed constant throughout the timestep to make this appear 
  # on the plot we have to do the following if there's fewer than 1000 timesteps
  fmintsp = 1000
  if integ_len < fmintsp:
    # work out small you need to divide each timestep to get 1000 timesteps
    div = ceil(fmintsp/integ_len)
    ftime = np.arange(0,integ_len,1.0/div) + y_0
    fluxes = ['emms','forc']
    for f in fluxes:
      tmp = []
      for j,v in enumerate(pts[f]):
        for i in range(0,int(div)):
          tmp.append(v)
      pts[f] = tmp
            
  else:
    ftime = time - 0.5
    
  if not infig:
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
  else:
    fig = infig
    ax1 = ax1in
    ax2 = ax2in
    ax3 = ax3in
    ax4 = ax4in

  ax1.plot(ftime,pts['emms'],color=colour['emms'],label=label,ls=linestyle)
  ax1.set_ylabel('Emissions (GtC)')
  if label is not None:
    ax1.legend(loc='best')
  ax2.plot(time,pts['conc'],color=colour['conc'],ls=linestyle)
  ax2.set_ylabel('CO$_2$ concentrations (ppm)')
  ax2.set_xlim(ax1.get_xlim())
  ax3.plot(ftime,pts['forc'],color=colour['forc'],ls=linestyle)
  ax3.set_ylabel('Other radiative forcing (W.m$^{-2}$)')
  ax3.set_xlabel('Time ({0})'.format(tuts))
  ax4.plot(time,pts['temp'],color=colour['temp'],ls=linestyle)
  ax4.set_ylabel('Temperature anomaly (K)')
  ax4.set_xlabel(ax3.get_xlabel())
  ax4.set_xlim(ax3.get_xlim())
  fig.tight_layout()

  return fig,ax1,ax2,ax3,ax4