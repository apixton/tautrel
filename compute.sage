attach "tautrel.sage"

ENABLE_DSAVE = True

def mod_type_string(moduli_type):
  if moduli_type == MODULI_ST:
    return "st"
  elif moduli_type == MODULI_CT:
    return "ct"
  elif moduli_type == MODULI_RT:
    return "rt"
  elif moduli_type == MODULI_SM:
    return "sm"
  return "xx"

def compute_betti(g,r,n=0,moduli_type=MODULI_ST):
  id_string = "betti-%s-%s-%s-%s" % (mod_type_string(moduli_type),g,r,n)
  start_time = time.time()
  start_memory = floor(get_memory_usage())
  #dsave("compute/%s-START",id_string)
  ans = betti(g,r,tuple([1 for i in range(n)]),moduli_type)
  end_time = time.time()
  end_memory = floor(get_memory_usage())
  dsave("compute/%s-ANSWER-%s-DUR-%s-MEM-%s",id_string,ans,floor(end_time-start_time),end_memory-start_memory)

def compute_goren(g,n=0,moduli_type=MODULI_ST):
  id_string = "goren-%s-%s-%s" % (mod_type_string(moduli_type),g,n)
  start_time = time.time()
  start_memory = floor(get_memory_usage())
  #dsave("compute/%s-START",id_string)
  d = dim_form(g,n,moduli_type)
  for r in range(0,floor(d/2)+1):
    ans = gorenstein(g,r,tuple([1 for i in range(n)]),moduli_type)
    end_time = time.time()
    end_memory = floor(get_memory_usage())
    dsave("compute/%s-ANSWER-%s-%s-DUR-%s-MEM-%s",id_string,r,ans,floor(end_time-start_time),end_memory-start_memory)

def compute_para_goren(g,n=0,moduli_type=MODULI_ST):
  id_string = "para_goren-%s-%s-%s" % (mod_type_string(moduli_type),g,n)
  start_time = time.time()
  start_memory = floor(get_memory_usage())
  #dsave("compute/%s-START",id_string)
  d = dim_form(g,n,moduli_type)
  for r in range(0,floor(d/2)+1):
    ans = para_gorenstein(g,r,tuple([1 for i in range(n)]),moduli_type)
    end_time = time.time()
    end_memory = floor(get_memory_usage())
    dsave("compute/%s-ANSWER-%s-%s-DUR-%s-MEM-%s",id_string,r,ans,floor(end_time-start_time),end_memory-start_memory)