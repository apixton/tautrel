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
  previous_rels = derived_rels(g,r,n,moduli_type)
  nrels = len(previous_rels)
  dlog('debug','choose_basic_rels(%s,%s,%s,%s): %s gens, %s oldrels, %s possible newrels',g,r,n,mod_type_string(moduli_type),sym_ngen,nrels,len(sym_possible_rels))
  D = {}
  for i in range(nrels):
    for x in previous_rels[i]:
      D[i,x[0]] = x[1]
  if nrels > 0:
    row_order,col_order = choose_orders_sparse(D,nrels,sym_ngen)
    previous_rank = compute_rank_sparse(D,row_order,col_order)
  else:
    previous_rank = 0
    row_order = []
    col_order = list(range(sym_ngen))
  answer = []
  for j in range(len(sym_possible_rels)):
    for x in sym_possible_rels[j]:
      D[nrels,x[0]] = x[1]
    row_order.append(nrels)
    nrels += 1
    if compute_rank_sparse(D,row_order,col_order) > previous_rank:
      answer.append(unsymmetrize_vec(sym_possible_rels[j],g,r,tuple(range(1,n+1)),moduli_type))
      previous_rank += 1
    if (j+1) % 5 == 0:
      dlog('debug','choose_basic_rels(%s,%s,%s,%s): checked %s newrels',g,r,n,mod_type_string(moduli_type),j)
  return answer

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
  rank = 0
  D = {}
  nrels = len(relations)
  if p > 0:
    KK = FiniteField(p)
  else:
    KK = QQ
  for i in range(nrels):
    for x in relations[i]:
      y = KK(x[1])
      if y != 0:
        D[i,x[0]] = y
  if nrels > 0:
    row_order,col_order = choose_orders_sparse(D,nrels,ngen)
    rank = compute_rank_sparse(D,row_order,col_order)
  return ngen - rank

def pullback_derived_rels(g,r,n=0,moduli_type=MODULI_ST):
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
            rel2 = insertion_pullback(rel2,g,r,n0+i,vec[i],moduli_type)
          answer.append(rel2)
    else:
      basic_rels = capply(choose_basic_rels,g,r,n0,MODULI_SMALL)
      k = r - dim_form(g,n0,moduli_type)
      for rel in basic_rels:
        for vec in subsequences(range(1,n0+k-1+1),k-1):
          for vec2 in subsequences(range(1,n+1),n-n0-k+1):
            rel2 = copy(rel)
            for i in range(k-1):
              rel2 = insertion_pullback(rel2,g,r,n0+i,vec[i],MODULI_SMALL)
            rel2 = insertion_pullback2(rel2,g,r,n0+k-1,vec2[0],moduli_type)
            for i in range(n-n0-k):
              rel2 = insertion_pullback(rel2,g,r,n0+k+i,vec2[i+1],moduli_type)
            answer.append(rel2)
  return answer

def interior_derived_rels(g,r,n=0,moduli_type=MODULI_ST):
  markings = tuple(range(1,n+1))
  answer = copy(capply(pullback_derived_rels,g,r,n,moduli_type))
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
                rel2 = psi_multiple(rel2,m+1,g,rcur,n,moduli_type)
                rcur += 1
            for m in sigma:
              rel2 = kappa_multiple(rel2,m,g,rcur,n,moduli_type)
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
