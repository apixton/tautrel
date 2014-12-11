attach "tautrel.sage"

def compute_betti(g,r,n=0,moduli_type=MODULI_ST):
  ans = log_func(betti,g,r,tuple([1 for i in range(n)]),moduli_type)

def compute_para_betti(g,r,n=0,moduli_type=MODULI_ST):
  ans = log_func(para_betti,g,r,tuple([1 for i in range(n)]),moduli_type)

def compute_goren(g,n=0,moduli_type=MODULI_ST):
  d = dim_form(g,n,moduli_type)
  for r in range(0,floor(d/2)+1):
    ans = log_func(gorenstein,g,r,tuple([1 for i in range(n)]),moduli_type)

def compute_para_goren(g,n=0,moduli_type=MODULI_ST):
  d = dim_form(g,n,moduli_type)
  for r in range(0,floor(d/2)+1):
    ans = log_func(para_gorenstein,g,r,tuple([1 for i in range(n)]),moduli_type)