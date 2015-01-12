class Graph:
  def __init__(self,M=None,genus_list=None):
    if M:
      self.M = copy(M)
    elif genus_list:
      self.M = matrix(R,len(genus_list)+1,1,[-1]+genus_list)
    else:
      self.M = matrix(R,1,1,-1)

  def num_vertices(self):
    return self.M.nrows() - 1

  def num_edges(self):
    return self.M.ncols() - 1

  def h1(self):
    return self.M.ncols()-self.M.nrows()+1

  def add_vertex(self,g):
    self.M = self.M.stack(matrix(1,self.M.ncols()))
    self.M[-1,0] = g

  def add_edge(self,i1,i2,marking=0):
    self.M = self.M.augment(matrix(self.M.nrows(),1))
    self.M[0,-1] = marking
    if i1 > 0:
      self.M[i1,-1] += 1
    if i2 > 0:
      self.M[i2,-1] += 1

  def del_vertex(self,i):
    self.M = self.M[:i].stack(self.M[(i+1):])

  def del_edge(self,i):
    if i == self.num_edges():
      self.M = self.M[0:,:i]
    else:
      self.M = self.M[0:,:i].augment(self.M[0:,(i+1):])

  def compute_degree_vec(self):
    self.degree_vec = [0 for i in range(1,self.M.nrows())]
    for i in range(1,self.M.nrows()):
      for j in range(1,self.M.ncols()):
        self.degree_vec[i-1] += self.M[i,j][0]

  def degree(self,i):
    return self.degree_vec[i-1]

  def split_vertex(self,i,row1,row2):
    self.M = self.M.stack(matrix(2,self.M.ncols(),row1+row2))
    self.add_edge(self.M.nrows()-2, self.M.nrows()-1)
    self.del_vertex(i)

  def set_target_parity(self):
    self.target_parity = 0
    for i in range(1,self.M.nrows()):
      local_parity = 1 + self.M[i,0][0]
      for j in range(1,self.M.ncols()):
        local_parity += self.M[i,j][1] + self.M[i,j][2]
      for j in range(1,self.M[i,0].degree()+1):
        local_parity += j*self.M[i,0][j]
      local_parity %= 2
      self.target_parity += (local_parity << (i-1))

  def replace_vertex_with_graph(self,i,G):
    nv = self.num_vertices()
    ne = self.num_edges()
    #i should have degree d, there should be no classes near i, and G should have markings 1,...,d and genus equal to the genus of i
    hedge_list = []
    for k in range(1,self.M.ncols()):
      for j in range(self.M[i,k]):
        hedge_list.append(k)
    self.del_vertex(i)
    for j in range(G.num_edges() - len(hedge_list)):
      self.add_edge(0,0)
    for j in range(G.num_vertices()):
      self.add_vertex(G.M[j+1,0])
    col = ne+1
    for k in range(1,G.M.ncols()):
      if G.M[0,k] > 0:
        mark = ZZ(G.M[0,k])
        for j in range(G.num_vertices()):
          if self.M[nv+j,hedge_list[mark-1]] == 0:
            self.M[nv+j,hedge_list[mark-1]] = G.M[j+1,k]
          elif G.M[j+1,k] != 0:
            a = self.M[nv+j,hedge_list[mark-1]][1]
            b = G.M[j+1,k][1]
            self.M[nv+j,hedge_list[mark-1]] = 2 + max(a,b)*X + min(a,b)*X^2
      else:
        for j in range(G.num_vertices()):
          self.M[nv+j,col] = G.M[j+1,k]
        col += 1

  def compute_invariant(self):
    nr,nc = self.M.nrows(),self.M.ncols()
    self.invariant = [[self.M[i,0], [], [], [[] for j in range(1,nr)]] for i in range(1,nr)]
    for k in range(1,nc):
      L = [i for i in range(1,nr) if self.M[i,k] != 0]
      if len(L) == 1:
        if self.M[0,k] != 0:
          self.invariant[L[0]-1][2].append((self.M[0,k],self.M[L[0],k]))
        else:
          self.invariant[L[0]-1][1].append(self.M[L[0],k])
      else:
        self.invariant[L[0]-1][3][L[1]-1].append((self.M[L[0],k],self.M[L[1],k]))
        self.invariant[L[1]-1][3][L[0]-1].append((self.M[L[1],k],self.M[L[0],k]))
    for i in range(1,nr):
      self.invariant[i-1][3] = [term for term in self.invariant[i-1][3] if len(term) > 0]
      for k in range(len(self.invariant[i-1][3])):
        self.invariant[i-1][3][k].sort()
        self.invariant[i-1][3][k] = tuple(self.invariant[i-1][3][k])
      self.invariant[i-1][3].sort()
      self.invariant[i-1][3] = tuple(self.invariant[i-1][3])
      self.invariant[i-1][2].sort()
      self.invariant[i-1][2] = tuple(self.invariant[i-1][2])
      self.invariant[i-1][1].sort()
      self.invariant[i-1][1] = tuple(self.invariant[i-1][1])
      self.invariant[i-1] = tuple(self.invariant[i-1])
    vertex_invariants = [[i,self.invariant[i-1]] for i in range(1,nr)]
    self.invariant.sort()
    self.invariant = tuple(self.invariant)
    vertex_invariants.sort(key=lambda x: x[1])
    self.vertex_groupings = []
    for i in range(nr-1):
      if i == 0 or vertex_invariants[i][1] != vertex_invariants[i-1][1]:
        self.vertex_groupings.append([])
      self.vertex_groupings[-1].append(vertex_invariants[i][0])

  def purify(self):
    for i in range(self.M.nrows()):
      for j in range(self.M.ncols()):
        self.M[i,j] = R(self.M[i,j][0])

  def contract(self,i,vlist,elist):
  # assumes graph is undecorated
    if self.M[0,i] != 0:
      print "ERROR: cannot contract a marking"
      return
    S = [row for row in range(1,self.M.nrows()) if self.M[row,i] != 0]
    if len(S) == 1:
      self.M[S[0],0] += 1
      self.del_edge(i)
      elist = elist[:(i-1)] + elist[i:]
    else:
      self.del_edge(i)
      elist = elist[:(i-1)] + elist[i:]
      self.add_vertex(0)
      self.M[-1] += self.M[S[0]]
      self.M[-1] += self.M[S[1]]
      self.del_vertex(S[1])
      self.del_vertex(S[0])
      vlist = vlist[:(S[0]-1)] + vlist[S[0]:(S[1]-1)] + vlist[S[1]:] + [vlist[S[0]-1] + vlist[S[1]-1]]
    return vlist,elist

def graph_isomorphic(G1,G2):
  if G1.invariant != G2.invariant:
    return False
  else:
    return isomorphic(G1.M,G2.M,G1.vertex_groupings,G2.vertex_groupings)

def isomorphic(M1,M2,group1,group2):
  nr,nc = M1.nrows(),M1.ncols()
  PermList = [Permutations(range(len(group))) for group in group1]
  for sigma_data in CartesianProduct(*PermList):
    sigma = [0 for i in range(nr-1)]
    for i in range(len(group1)):
      for j in range(len(group1[i])):
        sigma[group1[i][j]-1] = group2[i][sigma_data[i][j]]
    good = True
    for i in range(1,nr):
      ii = sigma[i-1]
      for j in range(1,i):
        jj = sigma[j-1]
        L1 = []
        for k in range(1,nc):
          if M1[i,k] != 0 and M1[j,k] != 0:
            L1.append([M1[i,k],M1[j,k]])
        L1.sort()
        L2 = []
        for k in range(1,nc):
          if M2[ii,k] != 0 and M2[jj,k] != 0:
            L2.append([M2[ii,k],M2[jj,k]])
        L2.sort()
        if L1 != L2:
          good = False
          break
      if good == False:
        break     
    if good:
      return True
  return False

def graph_count_automorphisms(G,vertex_orbits=False):
  return count_automorphisms(G.M,G.vertex_groupings,vertex_orbits)

def count_automorphisms(M,grouping,vertex_orbits=False):
  nr,nc = M.nrows(),M.ncols()
  count = 0
  PermList = [Permutations(range(len(group))) for group in grouping]
  if vertex_orbits:
    isom_list = []
  for sigma_data in CartesianProduct(*PermList):
    sigma = [0 for i in range(nr-1)]
    for i in range(len(grouping)):
      for j in range(len(grouping[i])):
        sigma[grouping[i][j]-1] = grouping[i][sigma_data[i][j]]
    good = True
    for i in range(1,nr):
      ii = sigma[i-1]
      for j in range(1,i):
        jj = sigma[j-1]
        L1 = []
        for k in range(1,nc):
          if M[i,k] != 0 and M[j,k] != 0:
            L1.append([M[i,k],M[j,k]])
        L1.sort()
        L2 = []
        for k in range(1,nc):
          if M[ii,k] != 0 and M[jj,k] != 0:
            L2.append([M[ii,k],M[jj,k]])
        L2.sort()
        if L1 != L2:
          good = False
          break
      if good == False:
        break     
    if good:
      count += 1
      if vertex_orbits:
        isom_list.append(sigma)

  if vertex_orbits:
    orbit_list = []
    vertices_used = []
    while len(vertices_used) < nr-1:
      i = [ii for ii in range(1,nr) if ii not in vertices_used][0]
      orbit = []
      for sigma in isom_list:
        if sigma[i-1] not in orbit:
          orbit.append(sigma[i-1])
          vertices_used.append(sigma[i-1])
      orbit.sort()
      orbit_list.append(orbit)
    return orbit_list

  for i in range(1,nr):
    for k in range(1,nc):
      if M[i,k][0] == 2 and M[i,k][1] == M[i,k][2]:
        count *= 2
    L = []
    for k in range(1,nc):
      if M[i,k] != 0:
        if sum(1 for j in range(1,nr) if M[j,k] != 0) == 1:
          L.append([M[0,k],M[i,k]])
    count *= aut(L)

    for j in range(1,i):
      L = []
      for k in range(1,nc):
        if M[i,k] != 0 and M[j,k] != 0:
          L.append([M[i,k],M[j,k]])
      count *= aut(L)
  return count

def graph_list_isomorphisms(G1,G2,only_one=False):
  if G1.invariant != G2.invariant:
    return []
  else:
    return list_isomorphisms(G1.M,G2.M,G1.vertex_groupings,G2.vertex_groupings,only_one)

def list_isomorphisms(M1,M2,group1,group2,only_one=False):
  # Warning: does not count loops!
  # If this is too slow, we can probably improve by caching a list of automorphisms and applying those to the first isom found.
  nr,nc = M1.nrows(),M2.ncols()
  PermList = [Permutations(range(len(group))) for group in group1]
  isom_list = []
  for sigma_data in CartesianProduct(*PermList):
    sigma = [0 for i in range(nr-1)]
    for i in range(len(group1)):
      for j in range(len(group1[i])):
        sigma[group1[i][j]-1] = group2[i][sigma_data[i][j]]
    good = True
    for i in range(1,nr):
      ii = sigma[i-1]
      for j in range(1,i):
        jj = sigma[j-1]
        L1 = []
        for k in range(1,nc):
          if M1[i,k] != 0 and M1[j,k] != 0:
            L1.append([M1[i,k],M1[j,k]])
        L1.sort()
        L2 = []
        for k in range(1,nc):
          if M2[ii,k] != 0 and M2[jj,k] != 0:
            L2.append([M2[ii,k],M2[jj,k]])
        L2.sort()
        if L1 != L2:
          good = False
          break
      if good == False:
        break     
    if good:
      cols1 = [[M1[i,j] for i in range(nr)] for j in range(1,nc)]
      cols2 = [[M2[0,j]] + [M2[sigma[i-1],j] for i in range(1,nr)] for j in range(1,nc)]
      edge_group1 = []
      edge_group2 = []
      used1 = []
      for j in range(1,nc):
        if j not in used1:
          edge_group1.append([])
          edge_group2.append([])
          for k in range(1,nc):
            if cols1[k-1] == cols1[j-1]:
              edge_group1[-1].append(k)
              used1.append(k)
            if cols2[k-1] == cols1[j-1]:
              edge_group2[-1].append(k)
      edge_PermList = [Permutations(range(len(edge_group))) for edge_group in edge_group1]
      for edge_sigma_data in CartesianProduct(*edge_PermList):
        edge_sigma = [0 for i in range(nc-1)]
        for i in range(len(edge_group1)):
          for j in range(len(edge_group1[i])):
            edge_sigma[edge_group1[i][j]-1] = edge_group2[i][edge_sigma_data[i][j]]
        isom_list.append([sigma,edge_sigma])
        if only_one:
          return isom_list
  return isom_list

def degenerate(G_list,moduli_type=MODULI_ST):
  mod_size = moduli_type + 1
  if moduli_type == MODULI_SMALL:
    mod_size = MODULI_SM + 1
  G_list_new = [[] for i in range(mod_size)]
  for which_type in range(mod_size):
    for G in G_list[which_type]:
      for i in range(1,G.num_vertices()+1):
        row = list(G.M[i])
        m = row[0] + sum(row)
        if m < 4:
          continue
        row1 = [0 for j in range(len(row))]
        while [2*x for x in row1] <= row:
          if row1[0] == 1 and moduli_type <= MODULI_RT:
            break
          if row1[0] + sum(row1) >= 2 and row1[0] + sum(row1) <= m-2:
            row2 = [row[j] - row1[j] for j in range(len(row))]
            G_copy = Graph(G.M)
            G_copy.split_vertex(i,row1,row2)
            new_type = which_type
            if new_type == MODULI_SM:
              new_type = MODULI_RT
            if new_type == MODULI_RT and row1[0] > 0:
              new_type = MODULI_CT
            G_list_new[new_type].append(G_copy)
          row1[-1] += 1
          for j in range(1,len(row)):
            if row1[-j] <= row[-j]:
              break
            row1[-j] = 0
            row1[-j-1] += 1
  for i in range(mod_size):
    G_list_new[i] = remove_isomorphic(G_list_new[i])
  return G_list_new

def dim_form(g,n,moduli_type=MODULI_ST):
  if moduli_type == MODULI_ST:
    return 3*g-3+n
  if moduli_type == MODULI_CT:
    return 2*g-3+n
  if moduli_type == MODULI_RT:
    if g > 0:
      return g-2+n
    else:
      return n-3
  if moduli_type == MODULI_SM:
    if n == 0:
      return g-2
    elif g >= 1:
      return g-1
    else:
      return 0
  if moduli_type == MODULI_SMALL:
    return 1000
  return 3*g-3+n

def decorate(G_list,r,moduli_type=MODULI_ST):
  mod_size = moduli_type + 1
  if moduli_type == MODULI_SMALL:
    mod_size = MODULI_SM + 1
  G_list_new = [[] for i in range(mod_size)]
  for which_type in range(mod_size):
    for G in G_list[which_type]:
      G_deco = [[] for i in range(mod_size)]
      G.compute_degree_vec()
      nr,nc = G.M.nrows(),G.M.ncols()
      two_list = []
      one_list = []
      for i in range(1,nr):
        for j in range(1,nc):
          if G.M[i,j] == 2:
            two_list.append([i,j])
          elif G.M[i,j] == 1:
            one_list.append([i,j])
      a = nr-1
      b = len(two_list)
      c = len(one_list)
      dims = [[dim_form(G.M[i+1,0][0], G.degree(i+1), mod_type) for i in range(a)] for mod_type in range(mod_size)]
      for vec in IntegerVectors(r,a+b+c):
        new_type = which_type
        if moduli_type > MODULI_SMALL:
          test_dims = vec[:a]
          for i in range(b):
            test_dims[two_list[i][0]-1] += vec[a+i]
          for i in range(c):
            test_dims[one_list[i][0]-1] += vec[a+b+i]
          for mod_type in range(which_type,mod_size):
            for i in range(a):
              if test_dims[i] > dims[mod_type][i]:
                new_type = mod_type + 1
                break
          if new_type > moduli_type:
            continue
        S_list = []
        for i in range(a):
          S_list.append(Partitions(vec[i]))
        for i in range(a,a+b):
          S_list.append([[vec[i]-j,j] for j in range(vec[i]/2 + 1)])
        S = CartesianProduct(*S_list)
        for vec2 in S:
          G_copy = Graph(G.M)
          for i in range(a):
            for j in vec2[i]:
              G_copy.M[i+1,0] += X^j
          for i in range(a,a+b):
            G_copy.M[two_list[i-a][0],two_list[i-a][1]] += vec2[i][0]*X + vec2[i][1]*X^2
          for i in range(c):
            G_copy.M[one_list[i][0],one_list[i][1]] += vec[i+a+b]*X
          G_deco[new_type].append(G_copy)
      for mod_type in range(mod_size):
        G_list_new[mod_type] += remove_isomorphic(G_deco[mod_type])
  return G_list_new

def remove_isomorphic(G_list):
  G_list_new = []
  inv_dict = {}
  count = 0
  for G1 in G_list:
    G1.compute_invariant()
    if not inv_dict.has_key(G1.invariant):
      inv_dict[G1.invariant] = []
    good = True
    for i in inv_dict[G1.invariant]:
      if graph_isomorphic(G1,G_list_new[i]):
        good = False
        break
    if good:
      G_list_new.append(G1)
      inv_dict[G1.invariant].append(count)
      count += 1
  return G_list_new

def num_strata(g,r,markings=(),moduli_type=MODULI_ST):
  return len(capply(all_strata,g,r,markings,moduli_type))

def num_pure_strata(g,r,markings=(),moduli_type=MODULI_ST):
  return len(capply(all_pure_strata,g,r,markings,moduli_type))

def single_stratum(num,g,r,markings=(),moduli_type=MODULI_ST):
  return capply(all_strata,g,r,markings,moduli_type)[num]

def single_pure_stratum(num,g,r,markings=(),moduli_type=MODULI_ST):
  return capply(all_pure_strata,g,r,markings,moduli_type)[num]

def autom_count(num,g,r,markings=(),moduli_type=MODULI_ST):
  return graph_count_automorphisms(single_stratum(num,g,r,markings,moduli_type))

def pure_strata_autom_count(num,g,r,markings=(),moduli_type=MODULI_ST):
  return graph_count_automorphisms(single_pure_stratum(num,g,r,markings,moduli_type))

def unpurify_map(g,r,markings=(),moduli_type=MODULI_ST):
  unpurify = {}
  pure_strata = [capply(all_pure_strata,g,r0,markings,moduli_type) for r0 in range(r+1)]
  impure_strata = capply(all_strata,g,r,markings,moduli_type)
  for i in range(len(impure_strata)):
    G = Graph(impure_strata[i].M)
    G.purify()
    r0 = G.num_edges() - len(markings)
    found = False
    for j in range(len(pure_strata[r0])):
      if G.M == pure_strata[r0][j].M:
        G_key = (r0, j)
        found = True
        break
    if not found:
      print "ERROR! Purification failed."
    if not unpurify.has_key(G_key):
      unpurify[G_key] = []
    unpurify[G_key].append(i)
  return unpurify

def all_strata(g,r,markings=(),moduli_type=MODULI_ST):
  dlog('debug','all_strata(%s,%s,%s,%s): begin',g,r,markings,mod_type_string(moduli_type))
  mod_size = moduli_type + 1
  if moduli_type == MODULI_SMALL:
    mod_size = MODULI_SM + 1
  big_list = [[] for i in range(mod_size)]
  for loops in range(g+1):
    if loops == 1 and moduli_type <= MODULI_CT:
      break
    if loops > r:
      break
    for edges in range(r-loops+1):
      if edges == 1 and moduli_type <= MODULI_SM:
        break
      G = Graph()
      G.add_vertex(g-loops)
      for k in range(loops):
        G.add_edge(1,1)
      for k in markings:
        G.add_edge(1,0,k)
      GGG = [[] for i in range(mod_size)]
      if loops == 0:
        if edges == 0:
          GGG[MODULI_SM] = [G]
        else:
          GGG[MODULI_RT] = [G]
      else:
        GGG[MODULI_ST] = [G]
      for k in range(edges):
        GGG = degenerate(GGG,moduli_type)
      GGG = decorate(GGG,r-loops-edges,moduli_type)
      for i in range(mod_size):
        big_list[i] += GGG[i]
  combined_list = []
  for i in range(mod_size):
    combined_list += big_list[i]
  for G in combined_list:
    G.compute_degree_vec()
    G.set_target_parity()
  dlog('debug','all_strata(%s,%s,%s,%s): %s gens',g,r,markings,mod_type_string(moduli_type),len(combined_list))
  return combined_list

def all_pure_strata(g,r,markings=(),moduli_type=MODULI_ST):
  big_list = [[] for i in range(moduli_type+1)]
  for loops in range(g+1):
    if loops == 1 and moduli_type <= MODULI_CT:
      break
    if loops > r:
      break
    for edges in range(r-loops,r-loops+1):
      if edges >= 1 and moduli_type <= MODULI_SM:
        break
      G = Graph()
      G.add_vertex(g-loops)
      for k in range(loops):
        G.add_edge(1,1)
      for k in markings:
        G.add_edge(1,0,k)
      G.compute_invariant()
      GGG = [[] for i in range(moduli_type+1)]
      if loops == 0:
        if edges == 0:
          GGG[MODULI_SM] = [G]
        else:
          GGG[MODULI_RT] = [G]
      else:
        GGG[MODULI_ST] = [G]
      for k in range(edges):
        GGG = degenerate(GGG,moduli_type)
      for i in range(moduli_type+1):
        big_list[i] += GGG[i]
  combined_list = []
  for i in range(moduli_type+1):
    combined_list += big_list[i]
  return combined_list

def strata_invariant_lookup(g,r,markings=(),moduli_type=MODULI_ST):
  inv_dict = {}
  L = capply(all_strata,g,r,markings,moduli_type)
  for i in range(len(L)):
    if not inv_dict.has_key(L[i].invariant):
      inv_dict[L[i].invariant] = []
    inv_dict[L[i].invariant].append(i)
  return inv_dict

def num_of_stratum(G,g,r,markings=(),moduli_type=MODULI_ST,special_top=False):
  if special_top:
    if g == TOP_g and r == TOP_r and markings == TOP_markings and moduli_type == TOP_moduli_type:
      return special_num_of_stratum(G,g,r,markings,moduli_type)
  G.compute_invariant()
  L = capply(all_strata,g,r,markings,moduli_type)
  x = capply(strata_invariant_lookup,g,r,markings,moduli_type)[G.invariant]
  if len(x) == 1:
    return x[0]
  for i in x:
    if graph_isomorphic(G,L[i]):
      return i
  print "ERROR"
  print (g,r,markings,moduli_type)
  print G.M

def special_num_of_stratum(G,g,r,markings=(),moduli_type=MODULI_ST):
  if markings == tuple([i+1 for i in range(len(TOP_new_markings))]):
    new_markings = TOP_new_markings
  else:
    print "ERROR"
    return "ERROR"
  GG = Graph(G.M)
  for i in range(1,GG.M.ncols()):
    if GG.M[0,i][0] > 0:
      GG.M[0,i] = R(new_markings[GG.M[0,i][0]-1])
  return num_of_stratum(GG,g,r,new_markings,moduli_type)