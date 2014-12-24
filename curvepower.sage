def all_strata_C(r,marktuple=()):
  markings = list(marktuple)
  n = len(markings)
  big_list = []
  for P in setparts(markings):
    kp_deg = r - n + len(P)
    if kp_deg < 0:
      continue
    block_sizes = []
    last = 0
    for i in range(1,len(P)):
      if P[i] != P[last]:
        block_sizes.append(i - last)
        last = i
    if last < len(P):
      block_sizes.append(len(P)-last)
    for deg_vec in IntegerVectors(kp_deg, 1+len(block_sizes)):
      S_list = [Partitions(deg_vec[0])]
      for i in range(len(block_sizes)):
        S_list.append(Partitions(deg_vec[i+1] + block_sizes[i], length=block_sizes[i]))
      for S in CartesianProduct(*S_list):
        stratum = [list(S[0]),[]]
        count = 0
        for i in range(len(block_sizes)):
          for j in range(block_sizes[i]):
            stratum[1].append([P[count],S[i+1][block_sizes[i]-1-j]-1])
            count += 1
        big_list.append(stratum)
  return big_list

def is_C_stratum(G):
  nr = G.M.nrows()
  nc = G.M.ncols()
  for i in range(1,nr):
    if G.M[i,0][0] != 0:
      continue
    deg = G.degree(i)
    dim = 0
    for j in range(1,G.M[i,0].degree()+1):
      dim += j*G.M[i,0][j]
    for j in range(1,nc):
      dim += G.M[i,j][1]
    if dim != deg - 3:
      return false
  return true
end

def convert_to_C(G):
  nr = G.M.nrows()
  nc = G.M.ncols()
  L = [i for i in range(1,nr) if G.M[i,0][0] != 0]
  if len(L) != 1:
    return 0
  v = L[0]
  E = [j for j in range(1,nc) if G.M[v,j][0] != 0]
  which_branch_vert = [0 for i in range(0,nr)]
  which_branch_edge = [0 for i in range(0,nc)]
  which_branch_vert[v] = -1
  which_branch_vert[0] = -1
  which_branch_edge[0] = -1
  for i in range(len(E)):
    which_branch_edge[E[i]] = i+1
  did_something = True
  while did_something:
    did_something = False
    for i in range(1,nr):
      for j in range(1,nc):
        if G.M[i,j][0] != 0:
          if which_branch_vert[i] > 0 and which_branch_edge[j] == 0:
            which_branch_edge[j] = which_branch_vert[i]
            did_something = True
          if which_branch_vert[i] == 0 and which_branch_edge[j] > 0:
            which_branch_vert[i] = which_branch_edge[j]
            did_something = True
  kappa_poly = G.M[v,0]
  kappa_list = []
  for i in range(1,kappa_poly.degree()+1):
    for j in range(kappa_poly[i]):
      kappa_list.append(i)
  kappa_list.reverse()
  stratum = [kappa_list, [[[G.M[0,k][0] for k in range(1,nc) if which_branch_edge[k] == i+1 and G.M[0,k][0] != 0],G.M[v,E[i]][1]] for i in range(len(E))]]
  for i in stratum[1]:
    i[0].sort()
  stratum[1].sort()

  factor = 1
  for i in range(1,nr):
    if i == v:
      continue
    kappa_poly = G.M[i,0]
    kappa_list = []
    for k in range(1,kappa_poly.degree()+1):
      for j in range(kappa_poly[k]):
        kappa_list.append(k+1)
    factor *= multinomial([G.M[i,j][1] for j in range(1,nc)] + kappa_list)
  return stratum,factor

def basic_C_rels(g,r,markings):
  generators = capply(all_strata,g,r,markings,MODULI_RT)
  gen_list = []
  for i in range(len(generators)):
    if is_C_stratum(generators[i]):
      gen_list.append(i)
  ngen = len(gen_list)
  C_generators = capply(all_strata_C,r,markings)
  C_ngen = len(C_generators)
  params = FZ_param_list(3*r-g-1,markings)
  C_conversion = []
  for k in range(ngen):
    strata,coeff = convert_to_C(generators[gen_list[k]])
    C_conversion.append([-1,coeff])
    for i in range(C_ngen):
      if C_generators[i] == strata:
        C_conversion[-1][0] = i
        break
  C_relations = []
  for FZ_param in params:
    relation = [0 for i in range(C_ngen)]
    for i in range(ngen):
      j = C_conversion[i][0]
      coeff = C_conversion[i][1]
      relation[j] += coeff*FZ_coeff(gen_list[i],FZ_param,g,r,markings,MODULI_RT)
    C_relations.append(relation)
  return C_relations

def list_all_FZ_C_general(g,r,markings,para=False):
  setpartlist = setparts(markings)
  n = len(markings)
  final_generators = capply(all_strata_C,r,markings)
  final_ngen = len(final_generators)
  final_relations = []

  # Sometimes do codim one differently.
  special_codim_one = False
  if n >= 5 and markings[-1] == 2:
    n1 = len([i for i in markings if i == 1])
    n2 = len([i for i in markings if i == 2])
    if n1 >= 2 and n2 >= 2:
      special_codim_one = True

  for codim in range(0,max(min(r,n),1)):
    if not (codim == 1 and special_codim_one):
      new_markings = []
      for i in range(n-2*codim):
        if markings[i] == markings[i+2*codim]:
          new_markings.append(markings[i])
      max_old_value = max([0] + new_markings)
      new_value = max_old_value + 1
      while len(new_markings) < n-codim:
        new_markings.append(new_value)
        new_value += 1
      new_markings = tuple(new_markings)
      C_generators = capply(all_strata_C,r-codim,new_markings)
      C_ngen = len(C_generators)
      if para:
        C_relations = para_basic_C_rels(g,r-codim,new_markings)
      else:
        C_relations = basic_C_rels(g,r-codim,new_markings)
    else:
      new_markings_list = []
      selector_dict = {}
      if n1 > 2:
        new_markings_list.append((n1-2)*[1] + n2*[2] + [3])
        selector_dict[1,1] = (0,[[1],[2],[1,1]])
      else:
        selector_dict[1,1] = (0,[[2],[2],[1,1]])
      if n2 != n1-1:
        new_markings_list.append((n1-1)*[1] + (n2-1)*[2] + [3])
        selector_dict[1,2] = (len(new_markings_list)-1,[[1],[2],[1,2]])
      else:
        selector_dict[1,2] = (0,[[2],[1],[1,2]])
      if n2 > 2 and n1 != n2 and n2 != n1+1:
        new_markings_list.append(n1*[1] + (n2-2)*[2] + [3])
        selector_dict[2,2] = (len(new_markings_list)-1,[[1],[2],[2,2]])
      elif n2 == 2:
        selector_dict[2,2] = (selector_dict[1,2][0],[[1],[1],[2,2]])
      elif n1 == n2:
        selector_dict[2,2] = (0,[[2],[1],[2,2]])
      elif n2 == n1+1:
        selector_dict[2,2] = (selector_dict[1,2][0],[[2],[1],[2,2]])
      for i in range(len(new_markings_list)):
        new_markings_list[i] = tuple(new_markings_list[i])
      C_generators_list = [capply(all_strata_C,r-codim,new_mark) for new_mark in new_markings_list]
      C_ngen_list = [len(C_gens) for C_gens in C_generators_list]
      if para:
        C_relations_list = [para_basic_C_rels(g,r-codim,new_mark) for new_mark in new_markings_list]
      else:
        C_relations_list = [basic_C_rels(g,r-codim,new_mark) for new_mark in new_markings_list]

    for setpart in setpartlist:
      if len(setpart) != n-codim:
        continue
      if codim == 1 and special_codim_one:
        for part in setpart:
          if len(part) == 2:
            i = selector_dict[tuple(part)][0]
            glue_key = selector_dict[tuple(part)][1]
            break
        new_markings = new_markings_list[i]
        C_generators = C_generators_list[i]
        C_ngen = C_ngen_list[i]
        C_relations = C_relations_list[i]
      else:
        setpart2 = [[i] for i in new_markings if i <= max_old_value]
        setpart3 = []
        count = 0
        for part in setpart:
          if count < len(setpart2) and part == setpart2[count]:
            count += 1
          else:
            setpart3.append(part)
        glue_key = [[i] for i in range(1,max_old_value+1)] + setpart3
      modified_C_relations = [[0 for i in range(final_ngen)] for j in range(len(C_relations))]
      for k in range(C_ngen):
        X = C_generators[k]
        strata = [X[0],[[[],i[1]] for i in X[1]]]
        for j in range(len(X[1])):
          for i in X[1][j][0]:
            strata[1][j][0] += glue_key[i-1]
          strata[1][j][0].sort()
        strata[1].sort()
        for i in range(final_ngen):
          if final_generators[i] == strata:
            which_gen = i
            break
        for j in range(len(C_relations)):
          modified_C_relations[j][which_gen] += C_relations[j][k]
      final_relations += modified_C_relations
  if len(final_relations) == 0:
    return [[0 for i in range(final_ngen)]]
  return final_relations

def list_all_FZ_C(g,r,n,para=False):
  markings = tuple(range(1,n+1))
  setpartlist = setparts(markings)
  final_generators = capply(all_strata_C,r,markings)
  final_ngen = len(final_generators)
  final_relations = []

  for codim in range(0,max(min(r,n),1)):
    new_markings = tuple(range(1,n-codim+1))
    C_generators = capply(all_strata_C,r-codim,new_markings)
    C_ngen = len(C_generators)
    if para:
      C_relations = para_basic_C_rels(g,r-codim,new_markings)
    else:
      C_relations = basic_C_rels(g,r-codim,new_markings)

    for setpart in setpartlist:
      if len(setpart) != n-codim:
        continue
      modified_C_relations = [[0 for i in range(final_ngen)] for j in range(len(C_relations))]
      for k in range(C_ngen):
        X = C_generators[k]
        strata = [X[0],[[[],i[1]] for i in X[1]]]
        for j in range(len(X[1])):
          for i in X[1][j][0]:
            strata[1][j][0] += setpart[i-1]
          strata[1][j][0].sort()
        strata[1].sort()
        for i in range(final_ngen):
          if final_generators[i] == strata:
            which_gen = i
            break
        for j in range(len(C_relations)):
          modified_C_relations[j][which_gen] += C_relations[j][k]
      final_relations += modified_C_relations
  if len(final_relations) == 0:
    return [[0 for i in range(final_ngen)]]
  return final_relations

def list_all_FZ_C_sym(g,r,n,para=False):
  markings = (1 for i in range(n))
  setpartlist = Partitions(n).list()
  final_generators = capply(all_strata_C,r,markings)
  final_ngen = len(final_generators)
  final_relations = []

  for codim in range(0,max(min(r,n),1)):
    num_ones = max(n-2*codim,0)
    new_markings = [1 for i in range(num_ones)] + range(2,n-codim-num_ones+2)
    new_markings = tuple(new_markings)
    C_generators = capply(all_strata_C,r-codim,new_markings)
    C_ngen = len(C_generators)
    if para:
      C_relations = para_basic_C_rels(g,r-codim,new_markings)
    else:
      C_relations = basic_C_rels(g,r-codim,new_markings)

    for setpart in setpartlist:
      if len(setpart) != n-codim:
        continue
      modified_C_relations = [[0 for i in range(final_ngen)] for j in range(len(C_relations))]
      for k in range(C_ngen):
        X = C_generators[k]
        strata = [X[0],[[[],i[1]] for i in X[1]]]
        for j in range(len(X[1])):
          for i in X[1][j][0]:
            if i == 1:
              strata[1][j][0].append(1)
            else:
              strata[1][j][0] += [1 for ii in range(setpart[i-2])]
          strata[1][j][0].sort()
        strata[1].sort()
        for i in range(final_ngen):
          if final_generators[i] == strata:
            which_gen = i
            break
        for j in range(len(C_relations)):
          modified_C_relations[j][which_gen] += C_relations[j][k]
      final_relations += modified_C_relations
  if len(final_relations) == 0:
    return [[0 for i in range(final_ngen)]]
  return final_relations

# todo: add mod p option in here
def betti_C(g,r,d=0):
  """
  This function returns the predicted rank of the codimension r grading
  of the tautological ring of the dth symmetric power of the universal
  curve over the moduli space of smooth genus g curves.

  g,r,d should be nonnegative integers, and g should be at least 2.

  EXAMPLE:

  - rank R^5(C^4_11)^(S_4) = 87
      sage: betti_C(11,5,4)
      87
  """
  if r > g+d-2:
    return 0
  L = list_all_FZ_C_sym(g,r,d)

  row_order,col_order = choose_orders(L)
  return (len(L[0]) - compute_rank2(L,row_order,col_order))

def betti_C_unsym(p,g,r,d=0):
  if r > g+d-2:
    return 0
  L = list_all_FZ_C(g,r,d)

  if p > 0:
    KK = FiniteField(p)
    for rel in L:
      for i in range(len(rel)):
        rel[i] = KK(rel[i])

  row_order,col_order = choose_orders(L)
  return (len(L[0]) - compute_rank2(L,row_order,col_order))

def betti_C_general(p,g,r,markings=()):
  if r > g+len(markings)-2:
    return 0
  L = list_all_FZ_C_general(g,r,markings)

  if p > 0:
    KK = FiniteField(p)
    for rel in L:
      for i in range(len(rel)):
        rel[i] = KK(rel[i])

  row_order,col_order = choose_orders(L)
  return (len(L[0]) - compute_rank2(L,row_order,col_order))

def para_betti_C(p,g,r,d=0):
  if r > g+d-2:
    return 0
  L = list_all_FZ_C_sym(g,r,d,True)

  if p > 0:
    KK = FiniteField(p)
    for rel in L:
      for i in range(len(rel)):
        rel[i] = KK(rel[i])

  row_order,col_order = choose_orders(L)
  return (len(L[0]) - compute_rank2(L,row_order,col_order))

def para_betti_C_unsym(p,g,r,d=0):
  if r > g+d-2:
    return 0
  L = list_all_FZ_C(g,r,d,True)

  if p > 0:
    KK = FiniteField(p)
    for rel in L:
      for i in range(len(rel)):
        rel[i] = KK(rel[i])

  row_order,col_order = choose_orders(L)
  return (len(L[0]) - compute_rank2(L,row_order,col_order))

def para_betti_C_general(p,g,r,markings=()):
  if r > g+len(markings)-2:
    return 0
  L = list_all_FZ_C_general(g,r,markings,True)

  if p > 0:
    KK = FiniteField(p)
    for rel in L:
      for i in range(len(rel)):
        rel[i] = KK(rel[i])

  row_order,col_order = choose_orders(L)
  return (len(L[0]) - compute_rank2(L,row_order,col_order))

@parallel
def para_FZ_subrels(gen_list,param_list,g,r,markings=(),moduli_type=MODULI_ST):
  rel_list = []
  for param in param_list:
    rel_list.append([])
    for i in gen_list:
      rel_list[-1].append(FZ_coeff(i,param,g,r,markings,moduli_type))
  return rel_list

def para_basic_C_rels(g,r,markings):
  generators = capply(all_strata,g,r,markings,MODULI_RT)
  gen_list = []
  for i in range(len(generators)):
    if is_C_stratum(generators[i]):
      gen_list.append(i)
  ngen = len(gen_list)
  C_generators = capply(all_strata_C,r,markings)
  C_ngen = len(C_generators)
  params = FZ_param_list(3*r-g-1,markings)
  C_conversion = []
  for k in range(ngen):
    strata,coeff = convert_to_C(generators[gen_list[k]])
    C_conversion.append([-1,coeff])
    for i in range(C_ngen):
      if C_generators[i] == strata:
        C_conversion[-1][0] = i
        break

  a1 = len(gen_list)
  k1 = min(a1,100)
  dlog('debug','parallelizing %s tasks in para_basic_C_rels(%s,%s,%s)',k1,g,r,markings)
  Slist = [gen_list[floor(i*a1/k1):floor((i+1)*a1/k1)] for i in range(k1)]
  input_list = []
  for S in Slist:
    input_list.append((S,params,g,r,markings,MODULI_RT))

  input_list = random_permutation(input_list)
  task_count = 0
  announce_every = 10

  C_relations = [[0 for i in range(C_ngen)] for j in range(len(params))]
  for FZ_rel_list in para_FZ_subrels(input_list):
    task_count += 1
    if task_count % announce_every == 0:
      dlog('debug','completed %s tasks in para_basic_C_rels(%s,%s,%s)',task_count,g,r,markings)
    for i in range(len(FZ_rel_list[0][0][0])):
      for num in range(len(gen_list)):
        if gen_list[num] == FZ_rel_list[0][0][0][i]:
          break
      j = C_conversion[num][0]
      coeff = C_conversion[num][1]
      for ii in range(len(params)):
        C_relations[ii][j] += coeff*FZ_rel_list[1][ii][i]
  return C_relations

######################### now gorenstein code

def multi2C(sigma):
  g = sum(sigma)+2
  term = factorial(2*g-3+len(sigma))
  term *= (2*g-1).multifactorial(2)
  term /= factorial(2*g-1)
  for i in sigma:
    term /= (2*i+1).multifactorial(2)
  return term

def kC_helper(sigma):
  total = 0
  for i in setparts_with_auts(sorted(sigma)):
    total += (-1)^(len(i[0])) * i[1] * multi2C([sum(j) for j in i[0]])
  return total

def kC(sigma):
  k0 = 2*sum(sigma) + 2
  sigma_reduced = [i for i in sigma if i != 0]
  zero_factor = k0^(len(sigma) - len(sigma_reduced))
  sigma_reduced.sort()
  return zero_factor*(-1)^(len(sigma_reduced))*capply(kC_helper,tuple(sigma_reduced))

def pairing_C(gen1,gen2):
  kappa_list = gen1[0] + gen2[0]
  diag_list = [[copy(r[0]),r[1]] for r in gen1[1] + gen2[1]]
  sign = 1
  while len(diag_list) > 0:
    c = diag_list[0][0][0]
    if c in diag_list[0][0][1:]:
      diag_list[0][1] += 1
      sign *= -1
      diag_list[0][0].remove(c)
      diag_list[0][0].remove(c)
    else:
      for ddd in diag_list[1:]:
        if c in ddd[0]:
          diag_list[0][1] += ddd[1]
          diag_list[0][0] += ddd[0]
          diag_list[0][0].remove(c)
          diag_list[0][0].remove(c)
          diag_list.remove(ddd)
          break
    if len(diag_list[0][0]) == 0:
      if diag_list[0][1] == 0:
        return 0
      kappa_list.append(diag_list[0][1]-1)
      diag_list = diag_list[1:]
  return sign*kC(kappa_list)

def pairing_C_sym(gen1,gen2,markings):
  gen1_copy = [gen1[0],[[copy(r[0]),r[1]] for r in gen1[1]]]
  gen2_copy = [gen2[0],[[copy(r[0]),r[1]] for r in gen2[1]]]
  block_sizes = []
  last = 0
  for i in range(1,len(markings)):
    if markings[i] != markings[last]:
      block_sizes.append([last,i])
      last = i
  if last < len(markings):
    block_sizes.append([last,len(markings)])
  for block in block_sizes:
    c = markings[block[0]]
    count = block[0]+1
    for i in range(len(gen1[1])):
      for j in range(len(gen1[1][i][0])):
        if gen1[1][i][0][j] == c:
          gen1_copy[1][i][0][j] = count
          count += 1
    count = block[0]+1
    for i in range(len(gen2[1])):
      for j in range(len(gen2[1][i][0])):
        if gen2[1][i][0][j] == c:
          gen2_copy[1][i][0][j] = count
          count += 1
  total = 0
  perm_list = [Permutations(range(block[0]+1,block[1]+1)) for block in block_sizes]
  gen2_copy2 = [gen2[0],[[copy(r[0]),r[1]] for r in gen2_copy[1]]]
  for perms in CartesianProduct(*perm_list):
    sigma = [0]
    for perm in perms:
      sigma += perm
    for i in range(len(gen2[1])):
      for j in range(len(gen2[1][i][0])):
        gen2_copy2[1][i][0][j] = sigma[gen2_copy[1][i][0][j]]
    total += pairing_C(gen1_copy,gen2_copy2)
  return total

def sym_orbit(gen,n):
  orbit = []
  for sigma in Permutations(range(1,n+1)):
    count = 0
    new_gen = [gen[0],[]]
    for i in range(len(gen[1])):
      new_gen[1].append([[],gen[1][i][1]])
      for j in range(len(gen[1][i][0])):
        new_gen[1][i][0].append(sigma[count])
        count += 1
      new_gen[1][i][0].sort()
    new_gen[1].sort()
    if new_gen not in orbit:
      orbit.append(new_gen)
  return orbit

def sym_orbit_one(gen,n):
  orbit = []
  count = 0
  new_gen = [gen[0],[]]
  for i in range(len(gen[1])):
    new_gen[1].append([[],gen[1][i][1]])
    for j in range(len(gen[1][i][0])):
      count += 1
      new_gen[1][i][0].append(count)
  return new_gen

def gor_betti(g,r,n):
  L1 = all_strata_C(r,tuple(range(1,n+1)))
  L2 = all_strata_C(g-2+n-r,tuple(range(1,n+1)))
  M = Matrix(QQ,len(L1),len(L2))
  for i in range(len(L1)):
    for j in range(len(L2)):
      M[i,j] = pairing_C(L1[i],L2[j])
  return M.rank()

def gor_betti_sym(g,r,n):
  L1 = all_strata_C(r,tuple([1 for i in range(n)]))
  L2 = all_strata_C(g-2+n-r,tuple([1 for i in range(n)]))
  M = Matrix(QQ,len(L1),len(L2))
  for i in range(len(L2)):
    L2[i] = sym_orbit_one(L2[i],n)
  for i in range(len(L1)):
    orbit = sym_orbit(L1[i],n)
    for gen in orbit:
      for j in range(len(L2)):
        M[i,j] += pairing_C(gen,L2[j])
  return M.rank()

def gorenstein_C(g,r,d=0):
  """
  This function returns the rank of the codimension r grading of the
  Gorenstein quotient of the tautological ring of the dth symmetric power
  of the universal curve over the moduli space of smooth genus g curves.

  g,r,d should be nonnegative integers, and g should be at least 2.

  EXAMPLE:

  - rank Gor^5(C^4_11)^(S_4) = 87
      sage: gorenstein_C(11,5,4)
      87
  """
  if r > g+d-2:
    return 0
  return gor_betti_sym(g,r,d)

def gorenstein_C_unsym(g,r,d=0):
  if r > g+d-2:
    return 0
  return gor_betti(g,r,d)

@parallel
def para_pairing_dict(S1,S2,g,r1,n):
  D = {}
  r2 = g-2+n-r1
  L1 = capply(all_strata_C,r1,tuple([1 for i in range(n)]))
  L2 = capply(all_strata_C,r2,tuple([1 for i in range(n)]))
  rep_dict = {}
  for i2 in S2:
    rep_dict[i2] = sym_orbit_one(L2[i2],n)
  for i1 in S1:
    orbit = sym_orbit(L1[i1],n)
    for i2 in S2:
      D[i1,i2] = 0
      for gen in orbit:
        D[i1,i2] += pairing_C(gen,rep_dict[i2])
  return D

def para_gorenstein_C(g,r1,n=0):
  r2 = g-2+n-r1
  # precompute stuff
  D = capply(all_strata_C,r1,tuple([1 for i in range(n)]))
  ngen1 = len(D)
  D = capply(all_strata_C,r2,tuple([1 for i in range(n)]))
  ngen2 = len(D)
  for sigma in Partitions(g-2):
    sigma_list = list(sigma)
    sigma_list.sort()
    D = capply(kC_helper,tuple(sigma_list))
  # probably should do more stuff here

  # should choose actual generators?
  S1 = list(range(ngen1))
  S2 = list(range(ngen2))

  # chop generator lists into blocks
  a1 = len(S1)
  a2 = len(S2)
  k1 = min(a1,30)
  k2 = min(a2,30)
  dlog('debug','parallelizing %s tasks in para_gorenstein_C(%s,%s,%s)',k1*k2,g,r1,n)
  S1list = [S1[floor(i*a1/k1):floor((i+1)*a1/k1)] for i in range(k1)]
  S2list = [S2[floor(i*a2/k2):floor((i+1)*a2/k2)] for i in range(k2)]
  input_list = []
  for T1 in S1list:
    for T2 in S2list:
      input_list.append((T1,T2,g,r1,n))
  input_list = random_permutation(input_list)
  result_list = list(para_pairing_dict(input_list))
  result_dict = {}
  for res in result_list:
    result_dict.update(res[1])

  M = [[0 for i2 in S2] for i1 in S1]
  for i1 in range(len(S1)):
    for i2 in range(len(S2)):
      M[i1][i2] = result_dict[S1[i1],S2[i2]]

  # compute rank; maybe should do this in a different way at some point.
  return matrix(M).rank()