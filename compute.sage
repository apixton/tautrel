attach "tautrel.sage"

def compute_betti(g,r,n=0,moduli_type=MODULI_ST):
  ans = log_func(betti,g,r,tuple([1 for i in range(n)]),moduli_type)

def compute_para_betti(p,g,r,n=0,moduli_type=MODULI_ST):
  ans = log_func(para_betti_prime,p,g,r,tuple([1 for i in range(n)]),moduli_type)

def compute_para_betti2(p,g,r,n=0,moduli_type=MODULI_ST):
  ans = log_func(para_betti_prime,p,g,r,tuple([i+1 for i in range(n)]),moduli_type)

def compute_goren(g,n=0,moduli_type=MODULI_ST):
  d = dim_form(g,n,moduli_type)
  for r in range(0,floor(d/2)+1):
    ans = log_func(gorenstein,g,r,tuple([1 for i in range(n)]),moduli_type)

def compute_para_goren(g,n=0,moduli_type=MODULI_ST):
  d = dim_form(g,n,moduli_type)
  for r in range(0,floor(d/2)+1):
    ans = log_func(para_gorenstein,g,r,tuple([1 for i in range(n)]),moduli_type)

def compute_para_goren2(g,n=0,moduli_type=MODULI_ST):
  d = dim_form(g,n,moduli_type)
  for r in range(0,floor(d/2)+1):
    ans = log_func(para_gorenstein,g,r,tuple([i+1 for i in range(n)]),moduli_type)

def compute_betti_C(g,r,n=0):
  ans = log_func(betti_C,g,r,n)

def compute_betti_C_unsym(g,r,n=0):
  ans = log_func(betti_C_unsym,g,r,n)

def compute_para_betti_C(p,g,r,n=0):
  ans = log_func(para_betti_C,p,g,r,n)

def compute_para_betti_C_unsym(p,g,r,n=0):
  ans = log_func(para_betti_C_unsym,p,g,r,n)

def compute_goren_C(g,n=0):
  d = g-2+n
  for r in range(0,floor(d/2)+1):
    ans = log_func(gorenstein_C,g,r,n)

def compute_goren_C_unsym(g,n=0):
  d = g-2+n
  for r in range(0,floor(d/2)+1):
    ans = log_func(gorenstein_C_unsym,g,r,n)

def compute_para_goren_C(g,n=0):
  d = g-2+n
  for r in range(0,floor(d/2)+1):
    ans = log_func(para_gorenstein_C,g,r,n)