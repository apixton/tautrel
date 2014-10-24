if 'initialized' not in globals():
  import time
  R.<X> = PolynomialRing(ZZ,1,order='lex')
  A_list = [factorial(6*n)/(factorial(3*n)*factorial(2*n)) for n in range(100)]
  B_list = [factorial(6*n+1)/((6*n-1)*factorial(3*n)*factorial(2*n)) for n in range(100)]
  MODULI_SMALL = -1
  MODULI_SM = 0
  MODULI_RT = 1
  MODULI_CT = 2
  MODULI_ST = 3
  cache_dict = {}
  ENABLE_DPRINT = False
  ENABLE_DSAVE = False
  initialized = True