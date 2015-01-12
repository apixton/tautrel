def num_new_rels(g,r,n=0,moduli_type=MODULI_ST):
  return len(capply(choose_basic_rels,g,r,n,moduli_type))

def choose_basic_rels(g,r,n=0,moduli_type=MODULI_ST):
  if 3*r < g+n+1:
    return []
  sym_ngen = num_strata(g,r,tuple([1 for i in range(n)]),moduli_type)
  if moduli_type == MODULI_SMALL and r > dim_form(g,n,MODULI_SM):
    sym_possible_rels = [[[i,1]] for i in range(sym_ngen)]
  else:
    sym_possible_rels = possibly_new_FZ(g,r,n,moduli_type)
  if len(sym_possible_rels) == 0:
    return []
  dlog('debug','choose_basic_rels(%s,%s,%s,%s): %s gens',g,r,n,mod_type_string(moduli_type),sym_ngen)
  previous_rels = derived_rels(g,r,n,moduli_type)
  nrels = len(previous_rels)
  dlog('debug','choose_basic_rels(%s,%s,%s,%s): %s gens, %s oldrels, %s possible newrels',g,r,n,mod_type_string(moduli_type),sym_ngen,nrels,len(sym_possible_rels))
  D = {}
  for i in range(nrels):
    for x in previous_rels[i]:
      D[i,x[0]] = x[1]
  if nrels > 0:
    previous_rank = compute_rank_sparse2(D,nrels,sym_ngen)
  else:
    previous_rank = 0
  dlog('debug','choose_basic_rels(%s,%s,%s,%s): initial rank is %s',g,r,n,mod_type_string(moduli_type),previous_rank)
  answer = []
  for j in range(len(sym_possible_rels)):
    for x in sym_possible_rels[j]:
      D[nrels,x[0]] = x[1]
    if not reduce_last_row(D,nrels,sym_ngen):
      new_rel = []
      for x in D.keys():
        if x[0] == nrels:
          new_rel.append([x[1],D[x]])
      answer.append(unsymmetrize_vec(new_rel,g,r,tuple(range(1,n+1)),moduli_type))
      previous_rank += 1
    nrels += 1
    if (j+1) % 5 == 0:
      dlog('debug','choose_basic_rels(%s,%s,%s,%s): checked %s rels',g,r,n,mod_type_string(moduli_type),j+1)
  dlog('debug','choose_basic_rels(%s,%s,%s,%s): found %s newrels',g,r,n,mod_type_string(moduli_type),len(answer))
  return answer

def assess_basic_rels(g,r,n=0,moduli_type=MODULI_ST):
  if 3*r < g+n+1:
    return []
  sym_ngen = num_strata(g,r,tuple([1 for i in range(n)]),moduli_type)
  if moduli_type == MODULI_SMALL and r > dim_form(g,n,MODULI_SM):
    sym_possible_rels = [[[i,1]] for i in range(sym_ngen)]
  else:
    sym_possible_rels = possibly_new_FZ(g,r,n,moduli_type)
  if len(sym_possible_rels) == 0:
    return []
  dlog('debug','assess_basic_rels(%s,%s,%s,%s): %s gens',g,r,n,mod_type_string(moduli_type),sym_ngen)
  previous_rels = derived_rels(g,r,n,moduli_type)
  nrels = len(previous_rels)
  dlog('debug','assess_basic_rels(%s,%s,%s,%s): %s gens, %s oldrels, %s possible newrels',g,r,n,mod_type_string(moduli_type),sym_ngen,nrels,len(sym_possible_rels))
  D = {}
  for i in range(nrels):
    for x in previous_rels[i]:
      D[i,x[0]] = x[1]
  if nrels > 0:
    previous_rank = compute_rank_sparse2(D,nrels,sym_ngen)
  else:
    previous_rank = 0
  dlog('debug','assess_basic_rels(%s,%s,%s,%s): initial rank is %s',g,r,n,mod_type_string(moduli_type),previous_rank)
  MM = Matrix(QQ,len(sym_possible_rels),sym_ngen)
  for j in range(len(sym_possible_rels)):
    DD = copy(D)
    for x in sym_possible_rels[j]:
      DD[nrels,x[0]] = x[1]
    _ = reduce_last_row(DD,nrels,sym_ngen)
    for i in range(sym_ngen):
      if DD.has_key((nrels,i)):
        MM[j,i] = DD[nrels,i]
  return list(MM.kernel().basis())

def recursive_betti(p,g,r,markings=(),moduli_type=MODULI_ST):
  n = len(markings)
  if r > dim_form(g,n,moduli_type):
    return 0
  ngen = num_strata(g,r,markings,moduli_type)
  dlog('debug','recursive_betti(%s,%s,%s,%s,%s): %s gens',p,g,r,markings,mod_type_string(moduli_type),ngen)
  relations = []
  partial_sym_map = capply(partial_symmetrize_map,g,r,markings,moduli_type)
  for rel in capply(choose_basic_rels,g,r,n,moduli_type):
    rel2 = []
    for x in rel:
      rel2.append([partial_sym_map[x[0]], x[1]])
    rel2 = simplify_sparse(rel2)
    relations.append(rel2)
  for rel in capply(interior_derived_rels,g,r,n,moduli_type):
    rel2 = []
    for x in rel:
      rel2.append([partial_sym_map[x[0]], x[1]])
    rel2 = simplify_sparse(rel2)
    relations.append(rel2)
  if moduli_type > MODULI_SM:
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
            rel_list = copy(capply(choose_basic_rels,g2,r-r0,d,moduli_type))
            rel_list += capply(interior_derived_rels,g2,r-r0,d,moduli_type)
            for rel0 in rel_list:
              relation = []
              for x in rel0:
                num = x[0]
                if which_gen_list[num] != -1:
                  relation.append([which_gen_list[num], x[1]])
              relation = simplify_sparse(relation)
              relations.append(relation)
  dlog('debug','recursive_betti(%s,%s,%s,%s,%s): %s rels',p,g,r,markings,mod_type_string(moduli_type),len(relations))
  relations = remove_duplicates2(relations)
  dlog('debug','recursive_betti(%s,%s,%s,%s,%s): %s distinct rels',p,g,r,markings,mod_type_string(moduli_type),len(relations))
  count = 0
  for rel in relations:
    count += len(rel)
  dlog('debug','recursive_betti(%s,%s,%s,%s,%s): %s nonzero entries',p,g,r,markings,mod_type_string(moduli_type),count)

  nrels = len(relations)
  if nrels == 0:
    return ngen
  if p > 0:
    KK = FiniteField(p)
    M = Matrix(KK,nrels,ngen,sparse=True)
    for i in range(nrels):
      for x in relations[i]:
        y = KK(x[1])
        if y != 0:
          M[i,x[0]] = y
    return ngen - M.rank()
  else:
    KK = QQ

  D = {}
  for i in range(nrels):
    for x in relations[i]:
      y = KK(x[1])
      if y != 0:
        D[i,x[0]] = y
  rank = compute_rank_sparse2(D,nrels,ngen)

  return ngen - rank

def fix_globals():
  global TOP_g
  global TOP_r
  global TOP_markings
  global TOP_new_markings
  global TOP_moduli_type
  TOP_g = None
  TOP_r = None
  TOP_markings = None
  TOP_new_markings = None
  TOP_moduli_type = None

def recursive_betti2(p,g,r,markings=(),moduli_type=MODULI_ST):
  global TOP_g
  global TOP_r
  global TOP_markings
  global TOP_new_markings
  global TOP_moduli_type
  TOP_g = g
  TOP_r = r
  n = len(markings)
  TOP_markings = tuple([i+1 for i in range(n)])
  TOP_new_markings = markings
  TOP_moduli_type = moduli_type
  n = len(markings)
  if r > dim_form(g,n,moduli_type):
    fix_globals()
    return 0
  ngen = num_strata(g,r,markings,moduli_type)
  dlog('debug','recursive_betti2(%s,%s,%s,%s,%s): %s gens',p,g,r,markings,mod_type_string(moduli_type),ngen)
  relations = []
  for FZ_param in FZ_param_list(3*r-g-1,markings):
    good = True
    for j in FZ_param[0]:
      if (j%3) != 1:
        good = False
    for j in FZ_param[1]:
      for jj in j[1]:
        if jj != 1:
          good = False
    if not good:
      continue
    relation = []
    for j in range(ngen):
      coeff = capply(FZ_coeff,j,FZ_param,g,r,markings,moduli_type)
      if coeff != 0:
        relation.append([j,coeff])
    relations.append(relation)
  dlog('debug','recursive_betti2(%s,%s,%s,%s,%s): %s possibly_new rels',p,g,r,markings,mod_type_string(moduli_type),len(relations))
  relations += interior_derived_rels(g,r,n,moduli_type,True)
  dlog('debug','recursive_betti2(%s,%s,%s,%s,%s): %s interior rels',p,g,r,markings,mod_type_string(moduli_type),len(relations))
  if moduli_type > MODULI_SM:
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
            rel_list = copy(capply(choose_basic_rels,g2,r-r0,d,moduli_type))
            rel_list += capply(interior_derived_rels,g2,r-r0,d,moduli_type)
            for rel0 in rel_list:
              relation = []
              for x in rel0:
                num = x[0]
                if which_gen_list[num] != -1:
                  relation.append([which_gen_list[num], x[1]])
              relation = simplify_sparse(relation)
              relations.append(relation)
  dlog('debug','recursive_betti2(%s,%s,%s,%s,%s): %s rels',p,g,r,markings,mod_type_string(moduli_type),len(relations))
  relations = remove_duplicates2(relations)
  dlog('debug','recursive_betti2(%s,%s,%s,%s,%s): %s distinct rels',p,g,r,markings,mod_type_string(moduli_type),len(relations))
  count = 0
  for rel in relations:
    count += len(rel)
  dlog('debug','recursive_betti2(%s,%s,%s,%s,%s): %s nonzero entries',p,g,r,markings,mod_type_string(moduli_type),count)

  nrels = len(relations)
  if nrels == 0:
    fix_globals()
    return ngen
  if p > 0:
    KK = FiniteField(p)
    M = Matrix(KK,nrels,ngen,sparse=True)
    for i in range(nrels):
      for x in relations[i]:
        y = KK(x[1])
        if y != 0:
          M[i,x[0]] = y
    fix_globals()
    return ngen - M.rank()
  else:
    KK = QQ

  D = {}
  for i in range(nrels):
    for x in relations[i]:
      y = KK(x[1])
      if y != 0:
        D[i,x[0]] = y
  rank = compute_rank_sparse2(D,nrels,ngen)

  fix_globals()
  return ngen - rank

def pullback_derived_rels(g,r,n=0,moduli_type=MODULI_ST,special_top=False):
  if r == 0:
    return []
  answer = []
  for n0 in range(n):
    if dim_form(g,n0,moduli_type) >= r:
      basic_rels = capply(choose_basic_rels,g,r,n0,moduli_type)
      for rel in basic_rels:
        for vec in subsequences(range(1,n+1),n-n0):
          rel2 = copy(rel)
          for i in range(n-n0):
            rel2 = insertion_pullback(rel2,g,r,n0+i,vec[i],moduli_type,special_top)
          answer.append(rel2)
    else:
      basic_rels = capply(choose_basic_rels,g,r,n0,MODULI_SMALL)
      k = r - dim_form(g,n0,moduli_type)
      for rel in basic_rels:
        for vec in subsequences(range(1,n0+k-1+1),k-1):
          for vec2 in subsequences(range(1,n+1),n-n0-k+1):
            rel2 = copy(rel)
            for i in range(k-1):
              rel2 = insertion_pullback(rel2,g,r,n0+i,vec[i],MODULI_SMALL,special_top)
            rel2 = insertion_pullback2(rel2,g,r,n0+k-1,vec2[0],moduli_type,special_top)
            for i in range(n-n0-k):
              rel2 = insertion_pullback(rel2,g,r,n0+k+i,vec2[i+1],moduli_type,special_top)
            answer.append(rel2)
  return answer

def interior_derived_rels(g,r,n=0,moduli_type=MODULI_ST,special_top=False):
  markings = tuple(range(1,n+1))
  answer = copy(capply(pullback_derived_rels,g,r,n,moduli_type,special_top))
  for r0 in range(r):
    pullback_rels = copy(capply(choose_basic_rels,g,r0,n,moduli_type))
    pullback_rels += capply(pullback_derived_rels,g,r0,n,moduli_type)
    for rel in pullback_rels:
      for i in range(r-r0+1):
        for sigma in Partitions(i):
          for tau in IntegerVectors(r-r0-i,n):
            rel2 = copy(rel)
            rcur = r0
            for m in range(n):
              for mm in range(tau[m]):
                rel2 = psi_multiple(rel2,m+1,g,rcur,n,moduli_type,special_top)
                rcur += 1
            for m in sigma:
              rel2 = kappa_multiple(rel2,m,g,rcur,n,moduli_type,special_top)
              rcur += m
            answer.append(rel2)
  return answer

def derived_rels(g,r,n=0,moduli_type=MODULI_ST):
  markings = tuple([1 for i in range(n)])
  generators = capply(all_strata,g,r,markings,moduli_type)
  ngen = len(generators)
  sym_map = capply(symmetrize_map,g,r,tuple(range(1,n+1)),moduli_type)
  answer = []
  for rel in capply(interior_derived_rels,g,r,n,moduli_type):
    rel2 = []
    for x in rel:
      rel2.append([sym_map[x[0]], x[1]])
    rel2 = simplify_sparse(rel2)
    answer.append(rel2)
  if moduli_type <= MODULI_SM:
    return answer
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
          rel_list = copy(capply(choose_basic_rels,g2,r-r0,d,moduli_type))
          rel_list += capply(interior_derived_rels,g2,r-r0,d,moduli_type)
          for rel0 in rel_list:
            relation = []
            for x in rel0:
              num = x[0]
              if which_gen_list[num] != -1:
                relation.append([which_gen_list[num], x[1]])
            relation = simplify_sparse(relation)
            answer.append(relation)
  answer = remove_duplicates2(answer)
  return answer

def list_num_new_rels(moduli_type=MODULI_ST):
  data = {}
  for key in cache_dict.keys():
    if key[0] == 'choose_basic_rels' and key[1][3] == moduli_type:
      data[tuple(list(key[1])[:3])] = len(cache_dict[key])
  key_list = data.keys()
  key_list.sort()
  for key in key_list:
    print "%s: %s" % (key, data[key])

def new_rel_locs(moduli_type=MODULI_ST):
  data = {}
  for key in cache_dict.keys():
    if key[0] == 'choose_basic_rels' and key[1][3] == moduli_type:
      data[tuple(list(key[1])[:3])] = len(cache_dict[key])
  key_list = data.keys()
  key_list.sort()
  for key in key_list:
    if data[key] > 0:
      print "%s: %s" % (key, data[key])
