if 'initialized' not in globals():
  import time
  R.<X> = PolynomialRing(ZZ,1,order='lex')
  MODULI_SMALL = -1
  MODULI_SM = 0
  MODULI_RT = 1
  MODULI_CT = 2
  MODULI_ST = 3
  cache_dict = {}
  initialized = True

  # these variables are dangerous, since they might cause capplied functions
  # to change their value if misused
  TOP_g = None
  TOP_r = None
  TOP_markings = None
  TOP_new_markings = None
  TOP_moduli_type = None