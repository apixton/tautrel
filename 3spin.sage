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

def dual_C_coeff(i,j,parity):
  total = 0
  k = parity % 2
  while (floor(k/3) <= i):
    if (k % 3) == 2:
      k += 2
      continue
    total += (-1)^(floor(k/3))*C_coeff(k,i)*C_coeff(-2-k,j)
    k += 2
  return total

def kappa_coeff(sigma,kappa_0,target_partition):
  total = 0
  num_ones = sum(1 for i in sigma if i == 1)
  for i in range(0,num_ones+1):
    for injection in Permutations(range(len(target_partition)),len(sigma)-i):
      term = binomial(num_ones,i)*binomial(kappa_0 + len(target_partition) + i-1, i)*factorial(i)
      for j in range(len(sigma)-i):
        term *= C_coeff(sigma[j+i],target_partition[injection[j]])
      for j in range(len(target_partition)):
        if j in injection:
          continue
        term *= C_coeff(0,target_partition[j])
      total += term
  total = (-1)^(len(target_partition)+len(sigma))*total/aut(list(target_partition))
  return total

def FZ_kappa_factor(num,sigma,g,r,markings=(),moduli_type=MODULI_ST):
  G = single_stratum(num,g,r,markings,moduli_type)
  L = []
  nv = G.num_vertices()
  for i in range(1,nv+1):
    L.append((2*G.M[i,0][0]+G.degree(i)-2,G.M[i,0]-G.M[i,0][0]))
  LL = []
  tau = []
  for i in range(nv):
    min = -1
    for j in range(nv):
      if (i == 0 or L[j] > LL[-1] or L[j] == LL[-1] and j > tau[-1]) and (min == -1 or L[j] < L[min]):
        min = j
    tau.append(min)
    LL.append(L[min])
  factor_dict = capply(FZ_kappa_factor2,tuple(LL),sigma)
  factor_vec = [0 for i in range(1 << nv)]
  for parity_key in factor_dict.keys():
    parity = 0
    for i in range(nv):
      if parity_key[i] == 1:
        parity += 1 << tau[i]
    factor_vec[parity] = factor_dict[parity_key]
  return factor_vec

def FZ_marking_factor(num,marking_vec,g,r,markings=(),moduli_type=MODULI_ST):
  G = single_stratum(num,g,r,markings,moduli_type)
  nv = G.num_vertices()
  ne = G.num_edges()
  num_parities = 2^nv
  PPP_list = []
  for marks in marking_vec:
    PPP_list.append(Permutations(marks[1]))
  PPP = CartesianProduct(*PPP_list)
  marking_factors = [0 for i in range(num_parities)]
  incident_vertices = []
  for mark_type in marking_vec:
    incident_vertices.append([])
    for k in range(1,ne+1):
      if G.M[0,k] == mark_type[0]:
        for i in range(1,nv+1):
          if G.M[i,k] != 0:
            incident_vertices[-1].append((i-1,G.M[i,k][1]))
            break
  for perms in PPP:
    parity = 0
    marking_factor = 1
    for marks_index in range(len(marking_vec)):
      for count in range(len(incident_vertices[marks_index])):
        marking_factor *= C_coeff(perms[marks_index][count],incident_vertices[marks_index][count][1])
        parity ^^= (perms[marks_index][count] % 2) << incident_vertices[marks_index][count][0]
    marking_factors[parity] += marking_factor
  return marking_factors

def FZ_kappa_factor2(L,sigma):
  nv = len(L)
  mmm = max((0,)+sigma)
  sigma_grouped = [0 for i in range(mmm)]
  for i in sigma:
    sigma_grouped[i-1] += 1
  S_list = []
  for i in sigma_grouped:
    S_list.append(IntegerVectors(i,nv))
  S = CartesianProduct(*S_list)
  kappa_factors = {}
  for parity in CartesianProduct(*[(0,1) for i in range(nv)]):
    kappa_factors[tuple(parity)] = 0
  for assignment in S:
    assigned_sigma = [[] for j in range(nv)]
    for i in range(mmm):
      for j in range(nv):
        for k in range(assignment[i][j]):
          assigned_sigma[j].append(i+1)
    sigma_auts = 1
    parity = [0 for i in range(nv)]
    kappa_factor = 1
    for j in range(nv):
      sigma_auts *= aut(assigned_sigma[j])
      parity[j] += sum(assigned_sigma[j])
      parity[j] %= 2
      kappa_factor *= capply(kappa_coeff,tuple(assigned_sigma[j]),L[j][0],poly_to_partition(L[j][1]))
    kappa_factors[tuple(parity)] += kappa_factor/sigma_auts
  return kappa_factors

def FZ_hedge_factor(num,g,r,markings=(),moduli_type=MODULI_ST):
  G = single_stratum(num,g,r,markings,moduli_type)
  nv = G.num_vertices()
  num_parities = 2^nv
  ne = G.num_edges()
  edge_list = []
  for k in range(1,ne+1):
    if G.M[0,k] == 0:
      edge_list.append([k])
      for i in range(1,nv+1):
        if G.M[i,k] != 0:
          edge_list[-1].append(i)
        if G.M[i,k][0] == 2:
          edge_list[-1].append(i)
  hedge_factors = [0 for i in range(num_parities)]
  for edge_parities in CartesianProduct(*[[0,1] for i in edge_list]):
    parity = 0
    for i in range(len(edge_list)):
      if edge_parities[i] == 1:
        parity ^^= 1 << (edge_list[i][1]-1)
        parity ^^= 1 << (edge_list[i][2]-1)
    hedge_factor = 1
    for i in range(len(edge_list)):
      if edge_list[i][1] == edge_list[i][2]:
        hedge_factor *= capply(dual_C_coeff,G.M[edge_list[i][1],edge_list[i][0]][1],G.M[edge_list[i][1],edge_list[i][0]][2],edge_parities[i] % 2)
      else:
        hedge_factor *= capply(dual_C_coeff,G.M[edge_list[i][1],edge_list[i][0]][1],G.M[edge_list[i][2],edge_list[i][0]][1],edge_parities[i] % 2)
    hedge_factors[parity] += hedge_factor
  return hedge_factors

def FZ_coeff(num,FZ_param,g,r,markings=(),moduli_type=MODULI_ST):
  sigma = FZ_param[0]
  marking_vec = FZ_param[1]
  G = single_stratum(num,g,r,markings,moduli_type)
  nv = G.num_vertices()
  graph_auts = capply(autom_count,num,g,r,markings,moduli_type)
  h1_factor = 2^G.h1()
  num_parities = 2^nv

  marking_factors = capply(FZ_marking_factor,num,marking_vec,g,r,markings,moduli_type)
  kappa_factors = capply(FZ_kappa_factor,num,sigma,g,r,markings,moduli_type)
  hedge_factors = capply(FZ_hedge_factor,num,g,r,markings,moduli_type)

  total = 0
  for i in range(num_parities):
    if marking_factors[i] == 0:
      continue
    for j in range(num_parities):
      total += marking_factors[i]*kappa_factors[j]*hedge_factors[i ^^ j ^^ G.target_parity]

  total /= h1_factor*graph_auts
  return total

def interior_FZ(g,r,markings=(),moduli_type=MODULI_ST):
  ngen = num_strata(g,r,markings,moduli_type)
  #print "%s generators" % ngen
  relations = []
  FZpl = FZ_param_list(3*r-g-1,markings)
  #print "%s codim 0 relations to compute" % len(FZpl)
  ccccc = 0
  for FZ_param in FZpl:
  #  if ccccc % 5 == 0:
  #    print "%s done" % ccccc
    ccccc += 1
    relation = [capply(FZ_coeff,i,FZ_param,g,r,markings,moduli_type) for i in range(ngen)]
    relations.append(relation)
  return relations

def possibly_new_FZ(g,r,n=0,moduli_type=MODULI_ST):
  m = 3*r-g-1-n
  if m < 0:
    return []
  dprint("Start FZ (%s,%s,%s,%s): %s",g,r,n,moduli_type,floor(get_memory_usage()))
  markings = tuple([1 for i in range(n)])
  ngen = num_strata(g,r,markings,moduli_type)
  relations = []
  for i in range(m+1):
    if m-i % 2 == 1:
      continue
    for sigma in Partitions(i):
      if len([j for j in sigma if j%3 != 1]) > 0:
        continue
      if n > 0:
        FZ_param = (tuple(sigma), ((1, markings),))
      else:
        FZ_param = (tuple(sigma), ())
      relation = []
      for j in range(ngen):
        coeff = capply(FZ_coeff,j,FZ_param,g,r,markings,moduli_type)
        if coeff != 0:
          relation.append([j,coeff])
      relations.append(relation)
  dprint("End FZ (%s,%s,%s,%s): %s",g,r,n,moduli_type,floor(get_memory_usage()))
  return relations
  

def boundary_FZ(g,r,markings=(),moduli_type=MODULI_ST):
  if moduli_type <= MODULI_SM:
    return []
  generators = capply(all_strata,g,r,markings,moduli_type)
  ngen = len(generators)
  #print "%s generators" % ngen
  relations = []
  old_count = 0
  for r0 in range(1,r):
    strata = capply(all_strata,g,r0,markings,moduli_type)
    for G in strata:
      vertex_orbits = graph_count_automorphisms(G,True)
      for i in [orbit[0] for orbit in vertex_orbits]:
        good = True
        for j in range(G.M.ncols()):
          if R(G.M[i,j][0]) != G.M[i,j]:
            good = False
            break
        if good:
          g2 = G.M[i,0][0]
          if 3*(r-r0) < g2 + 1:
            continue
          d = G.degree(i)
          if dim_form(g2,d,moduli_type) < r-r0:
            continue
          strata2 = capply(all_strata,g2,r-r0,tuple(range(1,d+1)),moduli_type)
          which_gen_list = [-1 for num in range(len(strata2))]
          for num in range(len(strata2)):
            G_copy = Graph(G.M)
            G_copy.replace_vertex_with_graph(i,strata2[num])
            which_gen_list[num] = num_of_stratum(G_copy,g,r,markings,moduli_type)
          rFZpl = reduced_FZ_param_list(G,i,g2,d,3*(r-r0)-g2-1)
          #print "Computing %s relations to insert at vertex %s into" % (len(rFZpl), i)
          #print G.M
          ccccc = 0
          for FZ_param in rFZpl:
            relation = [0 for k in range(ngen)]
            for num in range(len(strata2)):
              if which_gen_list[num] != -1:
                relation[which_gen_list[num]] += capply(FZ_coeff,num,FZ_param,g2,r-r0,tuple(range(1,d+1)),moduli_type)
            relations.append(relation)
            ccccc += 1
          #  if ccccc % 5 == 0:
          #    print "%s done" % ccccc
  return relations

def list_all_FZ(g,r,markings=(),moduli_type=MODULI_ST):
  relations = copy(capply(interior_FZ,g,r,markings,moduli_type))
  if moduli_type > MODULI_SM:
    relations += capply(boundary_FZ,g,r,markings,moduli_type)
  if len(relations) == 0:
    ngen = num_strata(g,r,markings,moduli_type)
    relations.append([0 for i in range(ngen)])
  return relations

def reduced_FZ_param_list(G,v,g,d,n):
  params = FZ_param_list(n,tuple(range(1,d+1)))
  graph_params = []
  M = matrix(R,2,d+1)
  M[0,0] = -1
  for i in range(1,d+1):
    M[0,i] = i
  for p in params:
    G_copy = Graph(G.M)
    M[1,0] = -g-1
    for j in p[0]:
      M[1,0] += X^j
    for i in range(1,d+1):
      M[1,i] = 1 + p[1][i-1][1][0]*X
    G_p = Graph(M)
    G_copy.replace_vertex_with_graph(v,G_p)
    graph_params.append([p,G_copy])
  params_reduced = []
  graphs_seen = []
  for x in graph_params:
    x[1].compute_invariant()
    good = True
    for GG in graphs_seen:
      if graph_isomorphic(x[1],GG):
        good = False
        break
    if good:
      graphs_seen.append(x[1])
      params_reduced.append(x[0])
  return params_reduced
    
def FZ_param_list(n,markings=()):
  if n < 0:
    return []
  final_list = []
  mmm = max((0,)+markings)
  markings_grouped = [0 for i in range(mmm)]
  for i in markings:
    markings_grouped[i-1] += 1
  markings_best = []
  for i in range(mmm):
    if markings_grouped[i] > 0:
      markings_best.append([i+1,markings_grouped[i]])
  for j in range(n/2 + 1):
    for n_vec in IntegerVectors(n-2*j,1+len(markings_best)):
      S_list = [[list(sigma) for sigma in Partitions(n_vec[0]).list() if sum(1 for l in sigma if (l % 3) == 2) == 0]]
      for i in range(len(markings_best)):
        S_list.append(Partitions(n_vec[i+1]+markings_best[i][1],length=markings_best[i][1]).list())
        S_list[-1] = [[k - 1 for k in sigma] for sigma in S_list[-1] if sum(1 for l in sigma if (l % 3) == 0) == 0]          
      for S in CartesianProduct(*S_list):
        final_list.append((tuple(S[0]),tuple([(markings_best[k][0],tuple(S[k+1])) for k in range(len(markings_best))])))
  return final_list

def FZ_matrix(g,r,markings=(),moduli_type=MODULI_ST):
  return matrix(list_all_FZ(g,r,markings,moduli_type))

def betti(g,r,marked_points=(),moduli_type=MODULI_ST):
  """
  This function returns the predicted rank of the codimension r grading
  of the tautological ring of the moduli space of stable genus g curves
  with marked points labeled by the multiset marked_points.

  g and r should be nonnegative integers and marked_points should be a
  tuple of positive integers.

  The parameter moduli_type determines which moduli space to use:
  - MODULI_ST: all stable curves (this is the default)
  - MODULI_CT: curves of compact type
  - MODULI_RT: curves with rational tails
  - MODULI_SM: smooth curves

  EXAMPLES:

  - rank R^3(bar{M}_2) = 1
      sage: betti(2,3)
      1

  - rank R^2(bar{M}_{2,3}) = 44
      sage: betti(2,2,(1,2,3))
      44

  - rank R^2(bar{M}_{2,3})^{S_3} = 20
      sage: betti(2,2,(1,1,1))
      20

  - rank R^2(bar{M}_{2,3})^{S_2} = 32 (S_2 interchanging markings 1 and 2)
      sage: betti(2,2,(1,1,2))
      32

  - rank R^2(M^c_4) = rank R^3(M^c_4) = 6
      sage: betti(4,2,(),MODULI_CT)
      6
      sage: betti(4,3,(),MODULI_CT)
      6

  - rank R^8(M^rt_{17,2})^(S_2) < R^9(M^rt_{17,2})^(S_2)
      sage: betti(17,8,(1,1),MODULI_RT)
      122
      sage: betti(17,9,(1,1),MODULI_RT)
      123

  - rank R^9(M_{20,1}) < rank R^10(M_{20,1})
      sage: betti(20,9,(1,),MODULI_SM)
      75
      sage: betti(20,10,(1,),MODULI_SM)
      76
  """
  L = list_all_FZ(g,r,marked_points,moduli_type)
  L.reverse()
  return (len(L[0]) - compute_rank(L))

def FZ_rels(g,r,markings=(),moduli_type=MODULI_ST):
  return span(FZ_matrix(g,r,markings,moduli_type)).basis()

##### some crude parallelization:

@parallel
def para_FZ_rels(g,r,g2,r0,d,markings=(),moduli_type=MODULI_ST):
  if r0 == 0:
    return interior_FZ(g,r,markings,moduli_type)
  if moduli_type <= MODULI_SM:
    return []
  if 3*(r-r0) < g2 + 1:
    return []
  if dim_form(g2,d,moduli_type) < r-r0:
    return []
  generators = capply(all_strata,g,r,markings,moduli_type)
  ngen = len(generators)
  relations = []
  old_count = 0
  strata = capply(all_strata,g,r0,markings,moduli_type)
  for G in strata:
    vertex_orbits = graph_count_automorphisms(G,True)
    for i in [orbit[0] for orbit in vertex_orbits]:
      if G.M[i,0][0] != g2:
        continue
      if G.degree(i) != d:
        continue
      good = True
      for j in range(G.M.ncols()):
        if R(G.M[i,j][0]) != G.M[i,j]:
          good = False
          break
      if good:
        strata2 = capply(all_strata,g2,r-r0,tuple(range(1,d+1)),moduli_type)
        which_gen_list = [-1 for num in range(len(strata2))]
        for num in range(len(strata2)):
          G_copy = Graph(G.M)
          G_copy.replace_vertex_with_graph(i,strata2[num])
          which_gen_list[num] = num_of_stratum(G_copy,g,r,markings,moduli_type)
        rFZpl = reduced_FZ_param_list(G,i,g2,d,3*(r-r0)-g2-1)
        for FZ_param in rFZpl:
          relation = [0 for k in range(ngen)]
          for num in range(len(strata2)):
            if which_gen_list[num] != -1:
              relation[which_gen_list[num]] += capply(FZ_coeff,num,FZ_param,g2,r-r0,tuple(range(1,d+1)),moduli_type)
          relations.append(relation)
  return relations

def para_betti_prime(p,g,r,markings=(),moduli_type=MODULI_ST):
  # precompute a few things?
  for r0 in range(r+1):
    D = capply(all_strata,g,r0,markings,moduli_type)
  # precompute vertex_orbits?

  input_list = []
  input_list.append((g,r,0,0,0,markings,moduli_type))
  degree_mult = 2
  if moduli_type < MODULI_ST:
    degree_mult = 1
  for r0 in range(1,r):
    for g2 in range(g,-1,-1):
      if 3*(r-r0) < g2 + 1:
        continue
      for d in range(1,degree_mult*r0+len(markings)+1):
        if dim_form(g2,d,moduli_type) < r-r0:
          continue
        input_list.append((g,r,g2,r0,d,markings,moduli_type))
  result_list = list(para_FZ_rels(input_list))
  relations = []
  for res in result_list:
    relations += res[1]
  if len(relations) == 0:
    return num_strata(g,r,markings,moduli_type)

  if p > 0:
    KK = FiniteField(p)
    return len(relations[0]) - matrix(KK,relations).rank()
  else:
    relations.reverse()
    return (len(relations[0]) - compute_rank(relations))
