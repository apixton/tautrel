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

# multiplication by kappa_i in the psi-pushforward basis
def dict_kappa_mult(D,i):
  DD = {}
  for sigma,coeff in D.iteritems():
    newsigma = list(sigma) + [i]
    newsigma.sort()
    newsigma = tuple(newsigma)
    if not DD.has_key(newsigma):
      DD[newsigma] = coeff
    else:
      DD[newsigma] += coeff
    for j in range(len(sigma)):
      if j > 0 and sigma[j] == sigma[j-1]:
        DD[newsigma] -= coeff
        continue
      newsigma = list(sigma)
      newsigma[j] += i
      newsigma.sort()
      newsigma = tuple(newsigma)
      if not DD.has_key(newsigma):
        DD[newsigma] = -coeff
      else:
        DD[newsigma] -= coeff
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

def dict_max_truncate(D,max):
  for sigma in D.keys():
    if base.sum(sigma) > max:
      del D[sigma]

def dict_parity_truncate(D,min,parity):
  for sigma in D.keys():
    s = base.sum(sigma)
    if s < min or (s % 2) != parity:
      del D[sigma]

def dict_simplify(D,g):
  for sigma in D.keys(): 
    count = 0
    for i in sigma:
      if i > 0:
        break
      count += 1
    if count > 0:
      newsigma = sigma[count:]
      val = D[sigma]*binomial(2*g-2+len(sigma)-1,count)*factorial(count)
      if not D.has_key(newsigma):
        D[newsigma] = val
      else:
        D[newsigma] += val
      del D[sigma]

def dict_eval(D,sigma,g,num_ones):
  ans = 0
  sorted_sigma = sorted(sigma)
  for i in range(num_ones+1):
    newsigma = tuple([0]*i + sorted_sigma)
    if D.has_key(newsigma):
      ans += D[newsigma]*binomial(2*g-2+len(newsigma)-1,i)*factorial(i)
  return ans

def tree_coeffs3(g,d,m,rel_matrix,rows_written,min_part,D,chosen_field):
  if (m % 2) == 0:
    vec = []
    for tau in Partitions(d):
      sorted_tau = tuple(sorted(tau))
      if D.has_key(sorted_tau):
        vec.append(D[sorted_tau])
      else:
        vec.append(chosen_field.zero())
    i = rows_written[0]
    for j in range(len(vec)):
      rel_matrix[i,j] = vec[j]
    rows_written[0] += 1
    if (rows_written[0] % 100) == 0:
      dlog('debug','%s rels computed in betti_mg(%s,%s)',rows_written[0],g,d)
  for i in range(min_part,floor(m/3)+1):
    dict_max_truncate(D,d-i)
    tree_coeffs3(g,d,m-3*i,rel_matrix,rows_written,i,dict_kappa_mult(D,i),chosen_field)

def tree_coeffs2(g,d,rel_matrix,rows_written,cur_sigma,D,chosen_field):
  s = sum(cur_sigma)
  maxm = 3*d-g-1-s
  if len(cur_sigma) > 0:
    minm = cur_sigma[-1]
  else:
    minm = 1
  for m in range(maxm,minm-1,-1):
    if (m % 3) != 1:
      continue
    L = [chosen_field(C_coeff(m,n)) for n in range(d+1)]
    tree_coeffs2(g,d,rel_matrix,rows_written,cur_sigma+[m],dict_mult(D,L,d),chosen_field)
  if maxm == 1:
    return
  if maxm < 5 and (maxm % 2) == 0:
    num_ones = sum(1 for i in cur_sigma if i == 1)
    vec = []
    for tau in Partitions(d):
      vec.append(dict_eval(D,list(tau),g,num_ones))
    i = rows_written[0]
    for j in range(len(vec)):
      rel_matrix[i,j] = vec[j]
    rows_written[0] += 1
    if (rows_written[0] % 100) == 0:
      dlog('debug','%s rels computed in betti_mg(%s,%s)',rows_written[0],g,d)
    return
  e = floor(maxm/3)
  mind = max(0,d-e)
  dict_parity_truncate(D,mind,(d-maxm)%2)
  dict_simplify(D,g)
  tree_coeffs3(g,d,maxm,rel_matrix,rows_written,1,D,chosen_field)

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
  if len(cur_sigma) > 0:
    minm = cur_sigma[-1]
  else:
    minm = 1
  for m in range(maxm,minm-1,-1):
    if (m % 3) == 2:
      continue
    L = [chosen_field(C_coeff(m,n)) for n in range(d+1)]
    tree_coeffs(g,d,rel_list,cur_sigma+[m],dict_mult(D,L,d),chosen_field)

def pnum(n,l):
  return Partitions(n,min_length=l).cardinality()

def sim_tree_coeffs(g,d,count_holder,cur_sigma):
  s = sum(cur_sigma)
  maxm = 3*d-g-1-s
  if len(cur_sigma) > 0:
    minm = cur_sigma[-1]
  else:
    minm = 1
  for m in range(maxm,minm-1,-1):
    if (m % 3) == 2:
      continue
    L = [QQ(C_coeff(m,n)) for n in range(d+1)]
    num_ones = sum(1 for i in cur_sigma if i == 1)
    for i in range(d+1):
      if L[i] != 0:
        for j in range(d-i+1):
          for k in range(num_ones + 1):
            count_holder[0] += capply(pnum,j,len(cur_sigma)-k)
    sim_tree_coeffs(g,d,count_holder,cur_sigma+[m])

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
  if p > 0:
    KK = FiniteField(p)
  else:
    KK = QQ
  nrows = mg_relcount(g,d)
  ncols = Partitions(d).cardinality()
  dlog('debug','computing %s rels in betti_mg(%s,%s)',mg_relcount(g,d),g,d)
  M = Matrix(KK,nrows,ncols)
  D = dict_exp_A(d,KK)
  tree_coeffs2(g,d,M,[0],[],D,KK)

  dlog('debug','computing rank in betti_mg(%s,%s,%s)',g,d,p)
  return ncols - M.rank()

def syz_mg(g,d,p=0):
  predicted_rank = Partitions(d).cardinality() - mg_relcount(g,d)
  actual_rank = betti_mg(g,d,p)
  return actual_rank - predicted_rank

def compute_betti_mg(g,r,p=0):
  ans = log_func(betti_mg,g,r,p)

def compute_syz_mg(g,r,p=0):
  ans = log_func(syz_mg,g,r,p)