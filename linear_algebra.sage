def compute_rank(L):
  count = 0
  for i in range(len(L)):
    S = [j for j in range(len(L[0])) if L[i][j] != 0]
    #print i,len(S)
    if len(S) == 0:
      continue
    count += 1
    j = S[0]
    T = [ii for ii in range(i+1,len(L)) if L[ii][j] != 0]
    for k in S[1:]:
      rat = L[i][k]/L[i][j]
      for ii in T:
        L[ii][k] -= rat*L[ii][j]
    for ii in range(i+1,len(L)):
      L[ii][j] = 0
  return count

def choose_orders_sparse(D,nrows,ncols):
  row_nums = [0 for i in range(nrows)]
  col_nums = [0 for j in range(ncols)]
  for key in D.keys():
    row_nums[key[0]] += 1
    col_nums[key[1]] += 1
  row_order = list(range(nrows))
  col_order = list(range(ncols))
  row_order.sort(key=lambda x: row_nums[x])
  col_order.sort(key=lambda x: col_nums[x])
  return row_order,col_order

def compute_rank_sparse(D,row_order,col_order):
  count = 0
  nrows = len(row_order)
  ncols = len(col_order)
  row_order_rank = [-1 for i in range(nrows)]
  col_order_rank = [-1 for i in range(ncols)]
  for i in range(nrows):
    row_order_rank[row_order[i]] = i
  for i in range(ncols):
    col_order_rank[col_order[i]] = i

  row_contents = [set() for i in range(nrows)]
  col_contents = [set() for i in range(ncols)]
  for x in D.keys():
    row_contents[x[0]].add(x[1])
    col_contents[x[1]].add(x[0])

  for i in row_order:
    S = []
    for j in row_contents[i]:
      S.append(j)
    if len(S) == 0:
      continue
    count += 1
    S.sort(key=lambda x: col_order_rank[x])
    j = S[0]
    T = []
    for ii in col_contents[j]:
      if row_order_rank[ii] > row_order_rank[i]:
        T.append(ii)
    for k in S[1:]:
      rat = D[i,k]/D[i,j]
      for ii in T:
        if not D.has_key((ii,k)):
          D[ii,k] = 0
          row_contents[ii].add(k)
          col_contents[k].add(ii)
        D[ii,k] -= rat*D[ii,j]
        if D[ii,k] == 0:
          D.pop((ii,k))
          row_contents[ii].remove(k)
          col_contents[k].remove(ii)
    for ii in T:
      D.pop((ii,j))
      row_contents[ii].remove(j)
      col_contents[j].remove(ii)
  return count

def compute_rank_sparse2(D,m,n):
  count = 0
  row_contents = [set() for i in range(m)]
  col_contents = [set() for i in range(n)]
  for x in D.keys():
    row_contents[x[0]].add(x[1])
    col_contents[x[1]].add(x[0])

  S = set([i for i in range(m) if len(row_contents[i]) > 0])

  while True:
    i = None
    ilen = n+1
    for ii in S:
      l = len(row_contents[ii])
      if l > 0 and l < ilen:
        i = ii
        ilen = l
    if i == None:
      return count
    count += 1
    j = None
    jlen = m+1
    for jj in row_contents[i]:
      l = len(col_contents[jj])
      if l < jlen:
        j = jj
        jlen = l
    for k in row_contents[i]:
      if k == j:
        continue
      rat = -D[i,k]/D[i,j]
      for ii in col_contents[j]:
        if ii == i:
          continue
        if not D.has_key((ii,k)):
          D[ii,k] = 0
          row_contents[ii].add(k)
          col_contents[k].add(ii)
        D[ii,k] += rat*D[ii,j]
        if D[ii,k] == 0:
          D.pop((ii,k))
          row_contents[ii].remove(k)
          col_contents[k].remove(ii)
    for ii in col_contents[j]:
      if ii == i:
        continue
      D.pop((ii,j))
      row_contents[ii].remove(j)
      if len(row_contents[ii]) == 0:
        S.remove(ii)
    for k in row_contents[i]:
      col_contents[k].remove(i)
    row_contents[i] = set()
    S.remove(i)
    col_contents[j] = set()

def choose_orders(L):
  rows = len(L)
  if rows == 0:
    return [],[]
  cols = len(L[0])
  row_nums = [0 for i in range(rows)]
  col_nums = [0 for j in range(cols)]
  for i in range(rows):
    for j in range(cols):
      if L[i][j] != 0:
        row_nums[i] += 1
        col_nums[j] += 1
  row_order = list(range(rows))
  col_order = list(range(cols))
  row_order.sort(key=lambda x: row_nums[x])
  col_order.sort(key=lambda x: col_nums[x])
  return row_order,col_order

def compute_rank2(L,row_order,col_order):
  count = 0
  for irow in range(len(row_order)):
    i = row_order[irow]
    S = [j for j in col_order if L[i][j] != 0]
    #print i,len(S)
    if len(S) == 0:
      continue
    count += 1
    j = S[0]
    T = [ii for ii in row_order[irow+1:] if L[ii][j] != 0]
    for k in S[1:]:
      rat = L[i][k]/L[i][j]
      for ii in T:
        L[ii][k] -= rat*L[ii][j]
    for ii in T:
      L[ii][j] = 0
  return count

def ranktest1(L,p=0):
  if p > 0:
    KK = FiniteField(p)
  else:
    KK = QQ
  m = len(L)
  n = 1 + max(max(x[0] for x in rel) for rel in L)
  D = {}
  for i in range(m):
    for x in L[i]:
      y = KK(x[1])
      if y != 0:
        D[i,x[0]] = y
  row_order,col_order = choose_orders_sparse(D,m,n)
  return compute_rank_sparse(D,row_order,col_order)

def ranktest2(L,p=0):
  if p > 0:
    KK = FiniteField(p)
  else:
    KK = QQ
  m = len(L)
  n = 1 + max(max(x[0] for x in rel) for rel in L)
  D = {}
  for i in range(m):
    for x in L[i]:
      y = KK(x[1])
      if y != 0:
        D[i,x[0]] = y
  return compute_rank_sparse2(D,m,n)

def ranktestsage(L,p=0):
  if p > 0:
    KK = FiniteField(p)
  else:
    KK = QQ
  m = len(L)
  n = 1 + max(max(x[0] for x in rel) for rel in L)
  M = Matrix(KK,m,n,sparse=True)
  for i in range(m):
    for x in L[i]:
      y = KK(x[1])
      if y != 0:
        M[i,x[0]] = y
  return M.rank()