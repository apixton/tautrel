import __builtin__ as base

attach "initialize.sage"
attach "meta.sage"
attach "util.sage"
attach "linear_algebra.sage"

def A_list(n):
  return factorial(6*n)/(factorial(3*n)*factorial(2*n))

def B_list(n):
  return factorial(6*n+1)/((6*n-1)*factorial(3*n)*factorial(2*n))

def C_coeff(m,term):
  n = term - floor(m/3)
  if n < 0:
    return 0
  if (m % 3) == 0:
    return capply(A_list,n)
  else:
    return capply(B_list,n)

def dict_mult(D,L,max_deg):
  DD = {}
  L_indices = [i for i in range(max_deg+1) if L[i] != 0]
  for sigma,coeff in D.iteritems():
    for i in L_indices:
      if base.sum(sigma)+i > max_deg:
        break
      newsigma = list(sigma) + [i]
      newsigma.sort()
      newsigma = tuple(newsigma)
      if not DD.has_key(newsigma):
        DD[newsigma] = L[i]*coeff
      else:
        DD[newsigma] += L[i]*coeff
  return DD

def dict_one():
  D = {}
  D[()] = 1
  return D

def dict_exp_A(max_deg, chosen_field):
  D = dict_one()
  D[()] = chosen_field(1)
  final_D = dict_one()
  final_D[()] = chosen_field(1)
  L = [-A_list(n) for n in range(max_deg+1)]
  L[0] = 0
  for n in range(max_deg+1):
    L[n] = chosen_field(L[n])
  for d in range(1,max_deg+1):
    D = dict_mult(D,L,max_deg)
    for key,val in D.iteritems():
      final_D[key] = val/factorial(d)
  return final_D

def dict_eval(D,sigma,g,num_ones):
  ans = 0
  sorted_sigma = sorted(sigma)
  for i in range(num_ones+1):
    newsigma = tuple([0]*i + sorted_sigma)
    if D.has_key(newsigma):
      ans += D[newsigma]*binomial(2*g-2+len(newsigma)-1,i)*factorial(i)
  return ans

def tree_coeffs(g,d,rel_list,cur_sigma,D,chosen_field):
  s = sum(cur_sigma)
  if ((s + g + d + 1) % 2) == 0:
    num_ones = sum(1 for i in cur_sigma if i == 1)
    vec = []
    for tau in Partitions(d):
      vec.append(dict_eval(D,list(tau),g,num_ones))
    rel_list.append(vec)
    if (len(rel_list) % 100) == 0:
      dlog('debug','%s rels computed in betti_mg(%s,%s)',len(rel_list),g,d)
  maxm = 3*d-g-1-s
  if len(cur_sigma) > 0 and maxm > cur_sigma[-1]:
    maxm = cur_sigma[-1]
  for m in range(1,maxm+1):
    if (m % 3) == 2:
      continue
    L = [chosen_field(C_coeff(m,n)) for n in range(d+1)]
    tree_coeffs(g,d,rel_list,cur_sigma+[m],dict_mult(D,L,d),chosen_field)

def mg_relcount(g,d):
  if 3*d-g-1 < 0:
    return 0
  count = 0
  for sigma in Partitions(3*d-g-1):
    good = True
    for part in sigma:
      if part > 2 and ((part % 3) == 2):
        good = False
        break
    if good:
      count += 1
  return count

def betti_mg(g,d,p=0):
  if d > g-2:
    return 0
  if 3*d < g+1:
    return Partitions(d).cardinality()
  dlog('debug','computing %s rels in betti_mg(%s,%s)',mg_relcount(g,d),g,d)
  rel_list = []
  if p > 0:
    KK = FiniteField(p)
  else:
    KK = QQ
  D = dict_exp_A(d,KK)
  tree_coeffs(g,d,rel_list,[],D,KK)

  dlog('debug','computing rank in betti_mg(%s,%s,%s)',g,d,p)
  M = Matrix(KK,rel_list)
  return M.ncols() - M.rank()
  #row_order,col_order = choose_orders(rel_list)
  #return (len(rel_list[0]) - compute_rank2(rel_list,row_order,col_order))

def syz_mg(g,d,p=0):
  predicted_rank = Partitions(d).cardinality() - mg_relcount(g,d)
  actual_rank = betti_mg(g,d,p)
  return actual_rank - predicted_rank

def compute_betti_mg(g,r,p=0):
  ans = log_func(betti_mg,g,r,p)

def compute_syz_mg(g,r,p=0):
  ans = log_func(syz_mg,g,r,p)