def kappa_conversion(sigma):
  answer = []
  for spart in setparts_with_auts(list(sigma)):
    coeff = spart[1]
    poly = R(0)
    for part in spart[0]:
      coeff *= factorial(len(part) - 1)
      poly += X^sum(part)
    answer.append([poly,coeff])
  return answer

def kappa_conversion_inverse(sigma):
  answer = []
  for spart in setparts_with_auts(list(sigma)):
    coeff = spart[1]*(-1)^(len(sigma)-len(spart[0]))
    poly = R(0)
    for part in spart[0]:
      poly += X^sum(part)
    answer.append([poly,coeff])
  return answer

def convert_to_monomial_basis(num,g,r,markings=(),moduli_type=MODULI_ST):
  answer = []
  G = single_stratum(num,g,r,markings,moduli_type)
  genus_vec = []
  kappa_vec = []
  for i in range(1,G.M.nrows()):
    genus_vec.append(G.M[i,0][0])
    kappa_vec.append([])
    for j in range(1,r+1):
      for k in range(G.M[i,0][j]):
        kappa_vec[-1].append(j)
    kappa_vec[-1] = capply(kappa_conversion,tuple(kappa_vec[-1]))
  for choice in CartesianProduct(*kappa_vec):
    coeff = 1
    GG = Graph(G.M)
    for i in range(1,G.M.nrows()):
      GG.M[i,0] = genus_vec[i-1] + choice[i-1][0]
      coeff *= choice[i-1][1]
    answer.append((num_of_stratum(GG,g,r,markings,moduli_type), coeff))
  return answer

def convert_to_pushforward_basis(num,g,r,markings=(),moduli_type=MODULI_ST):
  answer = []
  G = single_stratum(num,g,r,markings,moduli_type)
  genus_vec = []
  kappa_vec = []
  for i in range(1,G.M.nrows()):
    genus_vec.append(G.M[i,0][0])
    kappa_vec.append([])
    for j in range(1,r+1):
      for k in range(G.M[i,0][j]):
        kappa_vec[-1].append(j)
    kappa_vec[-1] = capply(kappa_conversion_inverse,tuple(kappa_vec[-1]))
  for choice in CartesianProduct(*kappa_vec):
    coeff = 1
    GG = Graph(G.M)
    for i in range(1,G.M.nrows()):
      GG.M[i,0] = genus_vec[i-1] + choice[i-1][0]
      coeff *= choice[i-1][1]
    answer.append((num_of_stratum(GG,g,r,markings,moduli_type), coeff))
  return answer

def convert_vector_to_monomial_basis(vec,g,r,markings=(),moduli_type=MODULI_ST):
  l = len(vec)
  vec2 = [0 for i in range(l)]
  for i in range(l):
    if vec[i] != 0:
      for x in capply(convert_to_monomial_basis,i,g,r,markings,moduli_type):
        vec2[x[0]] += x[1]*vec[i]
  return vec2

def convert_vector_to_pushforward_basis(vec,g,r,markings=(),moduli_type=MODULI_ST):
  l = len(vec)
  vec2 = [0 for i in range(l)]
  for i in range(l):
    if vec[i] != 0:
      for x in capply(convert_to_pushforward_basis,i,g,r,markings,moduli_type):
        vec2[x[0]] += x[1]*vec[i]
  return vec2

def kappa_multiple(vec,which_kappa,g,r,n=0,moduli_type=MODULI_ST,special_top=False):
  vec2 = []
  for x in vec:
    for y in capply(single_kappa_multiple,x[0],which_kappa,g,r,n,moduli_type,special_top):
      vec2.append([y[0], x[1]*y[1]])
  vec2 = simplify_sparse(vec2) 
  return vec2

def psi_multiple(vec,which_psi,g,r,n=0,moduli_type=MODULI_ST,special_top=False):
  vec2 = []
  for x in vec:
    for y in capply(single_psi_multiple,x[0],which_psi,g,r,n,moduli_type,special_top):
      vec2.append([y[0], x[1]*y[1]])
  vec2 = simplify_sparse(vec2)
  return vec2

def insertion_pullback(vec,g,r,n=0,new_mark=1,moduli_type=MODULI_ST,special_top=False):
  vec2 = []
  for x in vec:
    for y in capply(single_insertion_pullback,x[0],g,r,n,new_mark,moduli_type,special_top):
      vec2.append([y[0], x[1]*y[1]])
  vec2 = simplify_sparse(vec2)
  return vec2

def insertion_pullback2(vec,g,r,n=0,new_mark=1,moduli_type=MODULI_ST,special_top=False):
  vec2 = []
  for x in vec:
    for y in capply(single_insertion_pullback2,x[0],g,r,n,new_mark,moduli_type,special_top):
      vec2.append([y[0], x[1]*y[1]])
  vec2 = simplify_sparse(vec2)
  return vec2

def single_psi_multiple(num,which_psi,g,r,n=0,moduli_type=MODULI_ST,special_top=False):
  markings = tuple(range(1,n+1))
  G = single_stratum(num,g,r,markings,moduli_type)
  answer = []
  for j in range(1,G.M.ncols()):
    if G.M[0,j] == which_psi:
      good_j = j
      break
  for i in range(1,G.M.nrows()):
    if G.M[i,good_j] != 0:
      deg = 0
      dim_used = 0
      for j in range(1,r+1):
        dim_used += j*G.M[i,0][j]
      for j in range(1,G.M.ncols()):
        dim_used += G.M[i,j][1] + G.M[i,j][2]
        deg += G.M[i,j][0]
      if dim_used < dim_form(G.M[i,0][0],deg,moduli_type):
        GG = Graph(G.M)
        GG.M[i,good_j] += X
        answer.append((num_of_stratum(GG,g,r+1,markings,moduli_type,special_top),1))
      break
  return answer

def single_kappa_multiple(num,which_kappa,g,r,n=0,moduli_type=MODULI_ST,special_top=False):
  markings = tuple(range(1,n+1))
  G = single_stratum(num,g,r,markings,moduli_type)
  answer = []
  for i in range(1,G.M.nrows()):
    deg = 0
    dim_used = 0
    for j in range(1,r+1):
      dim_used += j*G.M[i,0][j]
    for j in range(1,G.M.ncols()):
      dim_used += G.M[i,j][1] + G.M[i,j][2]
      deg += G.M[i,j][0]
    if dim_used + which_kappa <= dim_form(G.M[i,0][0],deg,moduli_type):
      GG = Graph(G.M)
      GG.M[i,0] += X^which_kappa
      answer.append((num_of_stratum(GG,g,r+which_kappa,markings,moduli_type,special_top),1))
      for j in range(1,r+1):
        if G.M[i,0][j] > 0:
          GG = Graph(G.M)
          GG.M[i,0] += X^(j+which_kappa)
          GG.M[i,0] -= X^j
          answer.append((num_of_stratum(GG,g,r+which_kappa,markings,moduli_type,special_top),-G.M[i,0][j]))
  return answer

# Only doing this for markings=range(1,n+1) right now.
# Also, this function uses the monomial basis.
def single_kappa_psi_multiple(num,kappa_partition,psi_exps,g,r,n=0,moduli_type=MODULI_ST):
  markings = tuple(range(1,n+1))
  G = single_stratum(num,g,r,markings,moduli_type)
  GG = Graph(G.M)
  for j in range(1,G.M.ncols()):
    if GG.M[0,j] != 0:
      for i in range(1,G.M.nrows()):
        if GG.M[i,j] != 0:
          GG.M[i,j] += psi_exps[GG.M[0,j][0]-1]*X
          break
  rnew = r+sum(kappa_partition)+sum(psi_exps)
  answer = []
  kappa_options = [range(1,G.M.nrows()) for i in range(len(kappa_partition))]
  for kappa_distrib in CartesianProduct(*kappa_options):
    GGG = Graph(GG.M)
    for i in range(len(kappa_partition)):
      GGG.M[kappa_distrib[i],0] += X^(kappa_partition[i])
    is_bad = False
    for i in range(1,GGG.M.nrows()):
      deg = 0
      dim_used = 0
      for j in range(1,rnew+1):
        dim_used += j*GGG.M[i,0][j]
      for j in range(1,GGG.M.ncols()):
        dim_used += GGG.M[i,j][1] + GGG.M[i,j][2]
        deg += GGG.M[i,j][0]
      if dim_used > dim_form(GGG.M[i,0][0],deg,moduli_type):
        is_bad = True
        break
    if is_bad:
      continue
    answer.append((num_of_stratum(GGG,g,rnew,markings,moduli_type), 1))
  return answer

def single_insertion_pullback(num,g,r,n=0,new_mark=1,moduli_type=MODULI_ST,special_top=False):
  markings = tuple(range(1,n+1))
  new_markings = tuple(range(1,n+2))
  G = single_stratum(num,g,r,markings,moduli_type)
  answer = []
  for i in range(1,G.M.nrows()):
    GG = Graph(G.M)
    for j in range(1,G.M.ncols()):
      if GG.M[0,j][0] >= new_mark:
        GG.M[0,j] += 1
    GG.add_edge(i,0,new_mark)
    answer.append((num_of_stratum(GG,g,r,new_markings,moduli_type,special_top), 1))
    for j in range(1,r+1):
      for k in range(GG.M[i,0][j]):
        GGG = Graph(GG.M)
        GGG.M[i,0] -= X^j
        GGG.M[i,-1] += j*X
        answer.append((num_of_stratum(GGG,g,r,new_markings,moduli_type,special_top), -1))
    if moduli_type <= MODULI_SM:
      continue
    for j in range(1,G.M.ncols()):
      if G.M[i,j][0] == 1:
        if G.M[i,j][1] >= 1:
          x = G.M[i,j][1]
          GGG = Graph(GG.M)
          row1 = [GG.M[i,k] for k in range(GG.M.ncols())]
          row2 = [0 for k in range(GG.M.ncols())]
          row1[j] = 0
          row1[-1] = 0
          row2[j] = 1
          row2[-1] = 1
          GGG.split_vertex(i,row1,row2)
          GGG.M[-2,-1] += (x-1)*X
          answer.append((num_of_stratum(GGG,g,r,new_markings,moduli_type,special_top), -1))
      if G.M[i,j][0] == 2:
        if G.M[i,j][1] >= 1 or G.M[i,j][2] >= 1:
          x = G.M[i,j][1]
          y = G.M[i,j][2]
          row1 = [GG.M[i,k] for k in range(GG.M.ncols())]
          row2 = [0 for k in range(GG.M.ncols())]
          row1[j] = 0
          row1[-1] = 0
          row2[j] = 1
          row2[-1] = 1
          if y >= 1:
            row1[j] = 1 + x*X
            GGG = Graph(GG.M)
            GGG.split_vertex(i,row1,row2)
            GGG.M[-2,-1] += (y-1)*X
            answer.append((num_of_stratum(GGG,g,r,new_markings,moduli_type,special_top), -1))
          if x >= 1:
            row1[j] = 1 + y*X
            GGG = Graph(GG.M)
            GGG.split_vertex(i,row1,row2)
            GGG.M[-2,-1] += (x-1)*X
            answer.append((num_of_stratum(GGG,g,r,new_markings,moduli_type,special_top), -1))
  return answer

def single_insertion_pullback2(num,g,r,n=0,new_mark=1,moduli_type=MODULI_ST,special_top=False):
  markings = tuple(range(1,n+1))
  new_markings = tuple(range(1,n+2))
  G = single_stratum(num,g,r,markings,MODULI_SMALL)
  answer = []
  for i in range(1,G.M.nrows()):
    GG = Graph(G.M)
    for j in range(1,G.M.ncols()):
      if GG.M[0,j][0] >= new_mark:
        GG.M[0,j] += 1
    GG.add_edge(i,0,new_mark)
    answer.append((num_of_stratum(GG,g,r,new_markings,moduli_type,special_top), 1))
    for j in range(1,r+1):
      for k in range(GG.M[i,0][j]):
        GGG = Graph(GG.M)
        GGG.M[i,0] -= X^j
        GGG.M[i,-1] += j*X
        answer.append((num_of_stratum(GGG,g,r,new_markings,moduli_type,special_top), -1))
    if moduli_type <= MODULI_SM:
      continue
    for j in range(1,G.M.ncols()):
      if G.M[i,j][0] == 1:
        if G.M[i,j][1] >= 1:
          x = G.M[i,j][1]
          GGG = Graph(GG.M)
          row1 = [GG.M[i,k] for k in range(GG.M.ncols())]
          row2 = [0 for k in range(GG.M.ncols())]
          row1[j] = 0
          row1[-1] = 0
          row2[j] = 1
          row2[-1] = 1
          GGG.split_vertex(i,row1,row2)
          GGG.M[-2,-1] += (x-1)*X
          answer.append((num_of_stratum(GGG,g,r,new_markings,moduli_type,special_top), -1))
  return answer

def symmetrize_map(g,r,markings=(),moduli_type=MODULI_ST):
  markings2 = tuple([1 for i in markings])
  gens = capply(all_strata,g,r,markings,moduli_type)
  map = []
  for G in gens:
    GG = Graph(G.M)
    for i in range(1,GG.M.ncols()):
      if GG.M[0,i][0] > 0:
        GG.M[0,i] = R(1)
    map.append(num_of_stratum(GG,g,r,markings2,moduli_type))
  return map

def partial_symmetrize_map(g,r,markings=(),moduli_type=MODULI_ST):
  markings1 = tuple(range(1,len(markings)+1))
  gens = capply(all_strata,g,r,markings1,moduli_type)
  map = []
  for G in gens:
    GG = Graph(G.M)
    for i in range(1,GG.M.ncols()):
      if GG.M[0,i][0] > 0:
        GG.M[0,i] = R(markings[GG.M[0,i][0]-1])
    map.append(num_of_stratum(GG,g,r,markings,moduli_type))
  return map

def unsymmetrize_map(g,r,markings=(),moduli_type=MODULI_ST):
  markings2 = tuple([1 for i in markings])
  sym_map = capply(symmetrize_map,g,r,markings,moduli_type)
  map = [[] for i in range(num_strata(g,r,markings2,moduli_type))]
  for i in range(len(sym_map)):
    map[sym_map[i]].append(i)
  return map

def unsymmetrize_vec(vec,g,r,markings=(),moduli_type=MODULI_ST):
  unsym_map = capply(unsymmetrize_map,g,r,markings,moduli_type)
  vec2 = []
  for x in vec:
    aut = len(unsym_map[x[0]])
    for j in unsym_map[x[0]]:
      vec2.append([j,x[1]/aut])
  vec2 = simplify_sparse(vec2)
  return vec2