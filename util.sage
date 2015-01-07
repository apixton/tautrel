def aut(L):
  if len(L) == 0:
    return 1
  L.sort()
  total = 1
  n = 1
  last = L[0]
  for l in L[1:]:
    if l == last:
      n += 1
    else:
      n = 1
    total *= n
    last = l
  return total

def poly_to_partition(F):
  mmm = F.degree()
  target_partition = []
  for i in range(1,mmm+1):
    for j in range(F[i]):
      target_partition.append(i)
  return tuple(target_partition)

def setparts_recur(symlist,progress):
  if len(symlist) == 0:
    return [progress]
  l = []
  for i in Combinations(symlist[1:]).list():
    j = [symlist[0]]+i
    if len(progress) > 0 and j < progress[-1]:
      continue
    cur = 0
    new_symlist = []
    for k in range(len(symlist)):
      if cur < len(j) and symlist[k] == j[cur]:
        cur += 1
      else:
        new_symlist.append(symlist[k])
    l += setparts_recur(new_symlist,progress+[j])
  return l

def setparts_with_auts(symlist):
  l = setparts_recur(symlist,[])
  a = aut(symlist)
  ll = []
  for i in l:
    b = aut(i)
    for j in i:
      b *= aut(j)
    ll.append([i,a/b])
  return ll

def setparts(symlist):
  return setparts_recur(symlist,[])

def subsequences(seq,l):
  answer = []
  m = len(seq)
  for vec in IntegerVectors(l,m,max_part=1):
    answer.append([seq[i] for i in range(m) if vec[i] == 1])
  return answer

def submultisets_with_mults(symlist,l):
  if l < 0 or l > len(symlist):
    return []
  if len(symlist) == 0:
    return [[[],1,list(symlist)]]
  ans = []
  first_part = symlist[0]
  first_num = sum(1 for i in symlist if i == first_part)
  for i in range(first_num+1):
    trunc_list = symlist[first_num:]
    for x in submultisets_with_mults(trunc_list,l-i):
      ans.append([[first_part]*i+x[0],x[1]*binomial(first_num,i),[first_part]*(first_num-i)+x[2]])
  return ans

def injections_with_mults(L1,L2):
  if len(L1) > len(L2):
    return []
  if len(L1) == 0:
    return [[[],1,L2]]
  ans = []
  first_part = L1[0]
  first_num = 0
  for part in L1:
    if part != first_part:
      break
    first_num += 1
  LL1 = L1[first_num:]
  fact = factorial(first_num)
  for x in submultisets_with_mults(L2,first_num):
    partial_inj = [[first_part,a] for a in x[0]]
    for y in injections_with_mults(LL1,x[2]):
      ans.append([partial_inj + y[0],fact*x[1]*y[1],y[2]])
  return ans

def remove_duplicates(L):
  LL = []
  for elt in L:
    duplicate = False
    for elt2 in LL:
      if elt2 == elt:
        duplicate = True
        break
    if not duplicate:
      LL.append(elt)
  return LL

def remove_duplicates2(L):
  L.sort()
  LL = []
  for i in range(len(L)):
    if i == 0 or L[i] != L[i-1]:
      LL.append(L[i])
  return LL

def simplify_sparse(vec):
  vec.sort()
  vec2 = []
  last_index = None
  for x in vec:
    if x[0] == last_index:
      if vec2[-1][1] == -x[1]:
        vec2 = vec2[:-1]
        last_index = None
      else:
        vec2[-1][1] += x[1]
    elif x[1] != 0:
      vec2.append(x)
      last_index = x[0]
  return vec2

