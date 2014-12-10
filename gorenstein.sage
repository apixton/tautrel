def contraction_table(g,r,markings=(),moduli_type=MODULI_ST):
  contraction_dict = {}
  pure_strata = [capply(all_pure_strata,g,r0,markings,moduli_type) for r0 in range(r+1)]
  for r0 in range(r+1):
    for ii in range(len(pure_strata[r0])):
      G = pure_strata[r0][ii]
      S = [j for j in range(1,G.M.ncols()) if G.M[0,j] == 0]
      contractions = {}
      for edge_subset in CartesianProduct(*[[0,1] for i in range(r0)]):
        key = tuple(i for i in range(r0) if edge_subset[i] == 0)
        A = [S[i] for i in key]
        A.reverse()
        vlist = [[i] for i in range(1,G.M.nrows())]
        elist = [i for i in range(1,G.M.ncols())]
        Gcopy = Graph(G.M)
        for i in A:
          vlist,elist = Gcopy.contract(i,vlist,elist)
        Gcopy.compute_invariant()
        rnew = r0 - len(A)
        contraction_result = []
        for i in range(len(pure_strata[rnew])):
          L = graph_list_isomorphisms(pure_strata[rnew][i], Gcopy, True)
          if len(L) > 0:
            contraction_result.append((rnew,i))
            contraction_result.append(L[0])
            break
        contraction_result.append((vlist,elist))
        contractions[key] = contraction_result

      for edge_assignment in CartesianProduct(*[[0,1,2] for i in range(r0)]):
        if sum(1 for i in edge_assignment if i == 1) > r-r0:
          continue
        key1 = tuple(i for i in range(r0) if edge_assignment[i] == 0)
        B = [S[i] for i in range(r0) if edge_assignment[i] == 1]
        key2 = tuple(i for i in range(r0) if edge_assignment[i] == 2)
        if key1 > key2:
          continue
        contract1 = contractions[key1]
        contract2 = contractions[key2]
        dict_key = [contract1[0],contract2[0]]
        dict_entry = [contract1[1],contract2[1]]
        if dict_key[0] > dict_key[1]:
          dict_key.reverse()
          dict_entry.reverse()
          dict_entry = [(r0,ii),B,contract2[2],contract1[2]] + dict_entry
        else:
          dict_entry = [(r0,ii),B,contract1[2],contract2[2]] + dict_entry
        dict_key = tuple(dict_key)
        if not contraction_dict.has_key(dict_key):
          contraction_dict[dict_key] = []
        contraction_dict[dict_key].append(dict_entry)
        if dict_key[0] == dict_key[1]:
          contraction_dict[dict_key].append(dict_entry[:2] + [dict_entry[3],dict_entry[2],dict_entry[5],dict_entry[4]])
  return contraction_dict

def compute_pairing(r1,i1,r2,i2,g,rmax,markings=(),moduli_type=MODULI_ST):
  answer = 0
  pure_strata = [capply(all_pure_strata,g,r,markings,moduli_type) for r in range(rmax+1)]
  contraction_dict = capply(contraction_table,g,rmax,markings,moduli_type)
  G1 = single_stratum(i1,g,r1,markings,moduli_type)
  G2 = single_stratum(i2,g,r2,markings,moduli_type)
  G1copy = Graph(G1.M)
  G2copy = Graph(G2.M)
  G1copy.purify()
  G2copy.purify()
  pure_r1 = G1copy.num_edges() - len(markings)
  pure_r2 = G2copy.num_edges() - len(markings)
  found = False
  for i in range(len(pure_strata[pure_r1])):
    if G1copy.M == pure_strata[pure_r1][i].M:
      G1_key = (pure_r1, i)
      found = True
      break
  if not found:
    print "ERROR! Purification failed."
  found = False
  for i in range(len(pure_strata[pure_r2])):
    if G2copy.M == pure_strata[pure_r2][i].M:
      G2_key = (pure_r2, i)
      found = True
      break
  if not found:
    print "ERROR! Purification failed."
  if G1_key > G2_key:
    return compute_pairing(r2,i2,r1,i1,g,rmax,markings,moduli_type)

  if not contraction_dict.has_key((G1_key,G2_key)):
    return answer
  for L in contraction_dict[(G1_key,G2_key)]:
    H = pure_strata[L[0][0]][L[0][1]]
    Hloops = []
    if moduli_type > MODULI_CT:
      for i in range(1,H.M.nrows()):
        for j in range(1,H.M.ncols()):
          if H.M[i,j][0] == 2:
            Hloops.append((i,j))
    auts = capply(pure_strata_autom_count,L[0][1],g,L[0][0],markings,moduli_type)
    B = L[1]
    if len(B) == pure_r1 and len(B) == pure_r2:
      auts *= 2
    aut_cosets1 = capply(automorphism_cosets,i1,g,r1,markings,moduli_type)
    aut_cosets2 = capply(automorphism_cosets,i2,g,r2,markings,moduli_type)
    auts /= aut_cosets1[0]*aut_cosets2[0]
    for isom1 in aut_cosets1[1]:
      for isom2 in aut_cosets2[1]:
        Hcopy = Graph(H.M)
        vmap1 = [0 for i in range(G1.M.nrows())]
        for i in range(1,G1.M.nrows()):
          vmap1[i] = L[2][0][L[4][0][isom1[0][i-1]-1]-1]
        emap1 = [0 for i in range(G1.M.ncols())]
        for i in range(1,G1.M.ncols()):
          emap1[i] = L[2][1][L[4][1][isom1[1][i-1]-1]-1]
        vmap2 = [0 for i in range(G2.M.nrows())]
        for i in range(1,G2.M.nrows()):
          vmap2[i] = L[3][0][L[5][0][isom2[0][i-1]-1]-1]
        emap2 = [0 for i in range(G2.M.ncols())]
        for i in range(1,G2.M.ncols()):
          emap2[i] = L[3][1][L[5][1][isom2[1][i-1]-1]-1]

        psilooplist = []
        psiindexlist = []
        loop_factor = 1
        for i in range(1,G1.M.nrows()):
          for j in range(1,G1.M.ncols()):
            if G1.M[i,j][0] != 0:
              if G1.M[i,j][0] == 1:
                if G1.M[i,j][1] != 0:
                  jj = emap1[j]
                  for ii in vmap1[i]:
                    if H.M[ii,jj] != 0:
                      Hcopy.M[ii,jj] += G1.M[i,j][1]*X
                      break
              elif G1.M[i,j][1] == 0:
                loop_factor *= 2
              else:
                jj = emap1[j]
                psilooplist.append([[G1.M[i,j][1],G1.M[i,j][2]],[G1.M[i,j][2],G1.M[i,j][1]]])
                psiindexlist.append([jj])
                for ii in vmap1[i]:
                  for k in range(H.M[ii,jj][0]):
                    psiindexlist[-1].append(ii)
        for i in range(1,G2.M.nrows()):
          for j in range(1,G2.M.ncols()):
            if G2.M[i,j][0] != 0:
              if G2.M[i,j][0] == 1:
                if G2.M[i,j][1] != 0:
                  if G2.M[i,j][0] == 1:
                    jj = emap2[j]
                    for ii in vmap2[i]:
                      if H.M[ii,jj] != 0:
                        Hcopy.M[ii,jj] += G2.M[i,j][1]*X
                        break
              elif G2.M[i,j][1] == 0:
                loop_factor *= 2
              else:
                jj = emap2[j]
                psilooplist.append([[G2.M[i,j][1],G2.M[i,j][2]],[G2.M[i,j][2],G2.M[i,j][1]]])
                psiindexlist.append([jj])
                for ii in vmap2[i]:
                  for k in range(H.M[ii,jj][0]):
                    psiindexlist[-1].append(ii)

        Klocationlist = []
        Kindexlist = []
        for i in range(1,G1.M.nrows()):
          for r in range(1,rmax+1):
            for k in range(G1.M[i,0][r]):
              Klocationlist.append(vmap1[i])
              Kindexlist.append(r)
        for i in range(1,G2.M.nrows()):
          for r in range(1,rmax+1):
            for k in range(G2.M[i,0][r]):
              Klocationlist.append(vmap2[i])
              Kindexlist.append(r)

        psilist = []
        for j in B:
          S = [i for i in range(1,H.M.nrows()) if H.M[i,j][0] != 0]
          if len(S) == 2:
            psilist.append([[S[0],j],[S[1],j]])
          else:
            psilooplist.append([[0,1],[1,0]])
            psiindexlist.append([j,S[0],S[0]])

        for psiloopvals in CartesianProduct(*psilooplist):
          for Klocs in CartesianProduct(*Klocationlist):
            for psilocs in CartesianProduct(*psilist):
              Hcopycopy = Graph(Hcopy.M)
              for i in range(len(psiindexlist)):
                Hcopycopy.M[psiindexlist[i][1],psiindexlist[i][0]] += psiloopvals[i][0]*X
                if psiindexlist[i][1] == psiindexlist[i][2]:
                  Hcopycopy.M[psiindexlist[i][1],psiindexlist[i][0]] += psiloopvals[i][1]*X^2
                else:
                  Hcopycopy.M[psiindexlist[i][2],psiindexlist[i][0]] += psiloopvals[i][1]*X
              for i in range(len(Kindexlist)):
                Hcopycopy.M[Klocs[i],0] += X^Kindexlist[i]
              for i in psilocs:
                Hcopycopy.M[i[0],i[1]] += X
              for k in Hloops:
                if Hcopycopy.M[k][2] > Hcopycopy.M[k][1]:
                  Hcopycopy.M[k] += (Hcopycopy.M[k][2] - Hcopycopy.M[k][1])*(X - X^2)
              answer += (-1)^len(B)*loop_factor/auts*socle_eval(Hcopycopy,moduli_type)
  return answer

def multiply(r1,i1,r2,i2,g,rmax,markings=(),moduli_type=MODULI_ST):
  unpurify = capply(unpurify_map,g,r1+r2,markings,moduli_type)
  gens = capply(all_strata,g,r1+r2,markings,moduli_type)
  ngens = num_strata(g,r1+r2,markings,moduli_type)
  answer = [0 for i in range(ngens)]
  pure_strata = [capply(all_pure_strata,g,r,markings,moduli_type) for r in range(rmax+1)]
  contraction_dict = capply(contraction_table,g,rmax,markings,moduli_type)
  G1 = single_stratum(i1,g,r1,markings,moduli_type)
  G2 = single_stratum(i2,g,r2,markings,moduli_type)
  G1copy = Graph(G1.M)
  G2copy = Graph(G2.M)
  G1copy.purify()
  G2copy.purify()
  pure_r1 = G1copy.num_edges() - len(markings)
  pure_r2 = G2copy.num_edges() - len(markings)
  found = False
  for i in range(len(pure_strata[pure_r1])):
    if G1copy.M == pure_strata[pure_r1][i].M:
      G1_key = (pure_r1, i)
      found = True
      break
  if not found:
    print "ERROR! Purification failed."
  found = False
  for i in range(len(pure_strata[pure_r2])):
    if G2copy.M == pure_strata[pure_r2][i].M:
      G2_key = (pure_r2, i)
      found = True
      break
  if not found:
    print "ERROR! Purification failed."
  if G1_key > G2_key:
    return multiply(r2,i2,r1,i1,g,rmax,markings,moduli_type)

  if not contraction_dict.has_key((G1_key,G2_key)):
    return answer
  for L in contraction_dict[(G1_key,G2_key)]:
    H = pure_strata[L[0][0]][L[0][1]]
    Hloops = []
    if moduli_type > MODULI_CT:
      for i in range(1,H.M.nrows()):
        for j in range(1,H.M.ncols()):
          if H.M[i,j][0] == 2:
            Hloops.append((i,j))
    auts = capply(pure_strata_autom_count,L[0][1],g,L[0][0],markings,moduli_type)
    B = L[1]
    if len(B) == pure_r1 and len(B) == pure_r2:
      auts *= 2
    aut_cosets1 = capply(automorphism_cosets,i1,g,r1,markings,moduli_type)
    aut_cosets2 = capply(automorphism_cosets,i2,g,r2,markings,moduli_type)
    auts /= aut_cosets1[0]*aut_cosets2[0]
    for isom1 in aut_cosets1[1]:
      for isom2 in aut_cosets2[1]:
        Hcopy = Graph(H.M)
        vmap1 = [0 for i in range(G1.M.nrows())]
        for i in range(1,G1.M.nrows()):
          vmap1[i] = L[2][0][L[4][0][isom1[0][i-1]-1]-1]
        emap1 = [0 for i in range(G1.M.ncols())]
        for i in range(1,G1.M.ncols()):
          emap1[i] = L[2][1][L[4][1][isom1[1][i-1]-1]-1]
        vmap2 = [0 for i in range(G2.M.nrows())]
        for i in range(1,G2.M.nrows()):
          vmap2[i] = L[3][0][L[5][0][isom2[0][i-1]-1]-1]
        emap2 = [0 for i in range(G2.M.ncols())]
        for i in range(1,G2.M.ncols()):
          emap2[i] = L[3][1][L[5][1][isom2[1][i-1]-1]-1]

        psilooplist = []
        psiindexlist = []
        loop_factor = 1
        for i in range(1,G1.M.nrows()):
          for j in range(1,G1.M.ncols()):
            if G1.M[i,j][0] != 0:
              if G1.M[i,j][0] == 1:
                if G1.M[i,j][1] != 0:
                  jj = emap1[j]
                  for ii in vmap1[i]:
                    if H.M[ii,jj] != 0:
                      Hcopy.M[ii,jj] += G1.M[i,j][1]*X
                      break
              elif G1.M[i,j][1] == 0:
                loop_factor *= 2
              else:
                jj = emap1[j]
                psilooplist.append([[G1.M[i,j][1],G1.M[i,j][2]],[G1.M[i,j][2],G1.M[i,j][1]]])
                psiindexlist.append([jj])
                for ii in vmap1[i]:
                  for k in range(H.M[ii,jj][0]):
                    psiindexlist[-1].append(ii)
        for i in range(1,G2.M.nrows()):
          for j in range(1,G2.M.ncols()):
            if G2.M[i,j][0] != 0:
              if G2.M[i,j][0] == 1:
                if G2.M[i,j][1] != 0:
                  if G2.M[i,j][0] == 1:
                    jj = emap2[j]
                    for ii in vmap2[i]:
                      if H.M[ii,jj] != 0:
                        Hcopy.M[ii,jj] += G2.M[i,j][1]*X
                        break
              elif G2.M[i,j][1] == 0:
                loop_factor *= 2
              else:
                jj = emap2[j]
                psilooplist.append([[G2.M[i,j][1],G2.M[i,j][2]],[G2.M[i,j][2],G2.M[i,j][1]]])
                psiindexlist.append([jj])
                for ii in vmap2[i]:
                  for k in range(H.M[ii,jj][0]):
                    psiindexlist[-1].append(ii)

        Klocationlist = []
        Kindexlist = []
        for i in range(1,G1.M.nrows()):
          for r in range(1,rmax+1):
            for k in range(G1.M[i,0][r]):
              Klocationlist.append(vmap1[i])
              Kindexlist.append(r)
        for i in range(1,G2.M.nrows()):
          for r in range(1,rmax+1):
            for k in range(G2.M[i,0][r]):
              Klocationlist.append(vmap2[i])
              Kindexlist.append(r)

        psilist = []
        for j in B:
          S = [i for i in range(1,H.M.nrows()) if H.M[i,j][0] != 0]
          if len(S) == 2:
            psilist.append([[S[0],j],[S[1],j]])
          else:
            psilooplist.append([[0,1],[1,0]])
            psiindexlist.append([j,S[0],S[0]])

        for psiloopvals in CartesianProduct(*psilooplist):
          for Klocs in CartesianProduct(*Klocationlist):
            for psilocs in CartesianProduct(*psilist):
              Hcopycopy = Graph(Hcopy.M)
              for i in range(len(psiindexlist)):
                Hcopycopy.M[psiindexlist[i][1],psiindexlist[i][0]] += psiloopvals[i][0]*X
                if psiindexlist[i][1] == psiindexlist[i][2]:
                  Hcopycopy.M[psiindexlist[i][1],psiindexlist[i][0]] += psiloopvals[i][1]*X^2
                else:
                  Hcopycopy.M[psiindexlist[i][2],psiindexlist[i][0]] += psiloopvals[i][1]*X
              for i in range(len(Kindexlist)):
                Hcopycopy.M[Klocs[i],0] += X^Kindexlist[i]
              for i in psilocs:
                Hcopycopy.M[i[0],i[1]] += X
              for k in Hloops:
                if Hcopycopy.M[k][2] > Hcopycopy.M[k][1]:
                  Hcopycopy.M[k] += (Hcopycopy.M[k][2] - Hcopycopy.M[k][1])*(X - X^2)
              Hcopycopy.compute_invariant()
              for which_gen in unpurify[L[0]]:
                if graph_isomorphic(Hcopycopy,gens[which_gen]):
                  answer[which_gen] += (-1)^len(B)*loop_factor/auts
                  break
  return answer

def check_associativity(g,r1,r2,r3,markings=(),moduli_type=MODULI_ST):
  ngens1 = num_strata(g,r1,markings,moduli_type)
  ngens2 = num_strata(g,r2,markings,moduli_type)
  ngens3 = num_strata(g,r3,markings,moduli_type)
  print "%s*%s*%s = %s associators to compute:" % (ngens1, ngens2, ngens3, ngens1*ngens2*ngens3)
  count = 0
  for i1 in range(ngens1):
    for i2 in range(ngens2):
      for i3 in range(ngens3):
        a = capply(multiply,r1,i1,r2,i2,g,r1+r2+r3,markings,moduli_type)
        answer1 = vector([0 for i in range(num_strata(g,r1+r2+r3,markings,moduli_type))])
        for j in range(num_strata(g,r1+r2,markings,moduli_type)):
          if a[j] == 0:
            continue
          answer1 += a[j]*vector(capply(multiply,r1+r2,j,r3,i3,g,r1+r2+r3,markings,moduli_type))
        a = capply(multiply,r1,i1,r3,i3,g,r1+r2+r3,markings,moduli_type)
        answer2 = vector([0 for i in range(num_strata(g,r1+r2+r3,markings,moduli_type))])
        for j in range(num_strata(g,r1+r3,markings,moduli_type)):
          if a[j] == 0:
            continue
          answer2 += a[j]*vector(capply(multiply,r1+r3,j,r2,i2,g,r1+r2+r3,markings,moduli_type))
        if answer1 != answer2:
          print "Error: %s %s %s" % (i1,i2,i3)
        count += 1
        if count % 100 == 0:
          print "%s done" % count

def gorenstein_precompute(g,r1,markings=(),moduli_type=MODULI_ST):
  r3 = dim_form(g,len(markings),moduli_type)
  r2 = r3-r1
  D = capply(all_strata,g,r1,markings,moduli_type)
  D = capply(all_strata,g,r2,markings,moduli_type)
  D = capply(contraction_table,g,r3,markings,moduli_type)
  D = capply(unpurify_map,g,r3,markings,moduli_type)

def pairing_matrix(g,r1,markings=(),moduli_type=MODULI_ST):
  r3 = dim_form(g,len(markings),moduli_type)
  r2 = r3-r1
  ngens1 = num_strata(g,r1,markings,moduli_type)
  ngens2 = num_strata(g,r2,markings,moduli_type)
  ngens3 = num_strata(g,r3,markings,moduli_type)
  socle_evaluations = [socle_evaluation(i,g,markings,moduli_type) for i in range(ngens3)]
  pairings = [[0 for i2 in range(ngens2)] for i1 in range(ngens1)]
  if r1 == r2:
  #  print "%s*%s/2 = %s pairings to compute:" % (ngens1, ngens1+1, ngens1*(ngens1+1)/2)
    sym = True
  else:
  #  print "%s*%s = %s pairings to compute:" % (ngens1, ngens2, ngens1*ngens2)
    sym = False
  count = 0
  for i1 in range(ngens1):
    for i2 in range(ngens2):
      if sym and i1 > i2:
        pairings[i1][i2] = pairings[i2][i1]
        continue
      L = capply(multiply,r1,i1,r2,i2,g,r3,markings,moduli_type)
      pairings[i1][i2] = sum([L[k]*socle_evaluations[k] for k in range(ngens3)])
      count += 1
  #    if count % 100 == 0:
  #      print "%s done" % count
  return pairings

def pairing_submatrix(S1,S2,g,r1,markings=(),moduli_type=MODULI_ST):
  r3 = dim_form(g,len(markings),moduli_type)
  r2 = r3-r1
  ngens1 = num_strata(g,r1,markings,moduli_type)
  ngens2 = num_strata(g,r2,markings,moduli_type)
  ngens3 = num_strata(g,r3,markings,moduli_type)
  #print "Computing socle evaluation"
  socle_evaluations = [socle_evaluation(i,g,markings,moduli_type) for i in range(ngens3)]
  pairings = [[0 for i2 in S2] for i1 in S1]
  if r1 == r2 and S1 == S2:
  #  print "%s*%s/2 = %s pairings to compute:" % (len(S1), len(S1)+1, len(S1)*(len(S1)+1)/2)
    sym = True
  else:
  #  print "%s*%s = %s pairings to compute:" % (len(S1), len(S2), len(S1)*len(S2))
    sym = False
  count = 0
  for i1 in range(len(S1)):
    for i2 in range(len(S2)):
      if sym and i1 > i2:
        pairings[i1][i2] = pairings[i2][i1]
        continue
      L = capply(multiply,r1,S1[i1],r2,S2[i2],g,r3,markings,moduli_type)
      pairings[i1][i2] = sum([L[k]*socle_evaluations[k] for k in range(ngens3)])
      count += 1
  #    if count % 100 == 0:
  #      print "%s done" % count
  return pairings

def socle_eval(G,moduli_type=MODULI_ST):
  answer = QQ(1)
  for i in range(1,G.M.nrows()):
    g0 = G.M[i,0][0]
    psilist = []
    for j in range(1,G.M.ncols()):
      if G.M[i,j][0] > 0:
        psilist.append(G.M[i,j][1])
        if G.M[i,j][0] == 2:
          psilist.append(G.M[i,j][2])
    n0 = len(psilist)
    dim0 = dim_form(g0,n0,moduli_type)
    kappalist = []
    for j in range(1,dim0+1):
      for k in range(G.M[i,0][j]):
        kappalist.append(j)
    if sum(psilist)+sum(kappalist) != dim0:
      return QQ(0)
    answer *= capply(socle_formula,g0,tuple(psilist),tuple(kappalist),moduli_type)
  return answer

def socle_evaluation(num,g,markings=(),moduli_type=MODULI_ST):
  answer = 1
  G = single_stratum(num,g,dim_form(g,len(markings),moduli_type),markings,moduli_type)
  for i in range(1,G.M.nrows()):
    g0 = G.M[i,0][0]
    psilist = []
    for j in range(1,G.M.ncols()):
      if G.M[i,j][0] > 0:
        psilist.append(G.M[i,j][1])
        if G.M[i,j][0] == 2:
          psilist.append(G.M[i,j][2])
    n0 = len(psilist)
    dim0 = dim_form(g0,n0,moduli_type)
    kappalist = []
    for j in range(1,dim0+1):
      for k in range(G.M[i,0][j]):
        kappalist.append(j)
    if sum(psilist)+sum(kappalist) != dim0:
      print "ERROR: wrong dim"
      return
    answer *= capply(socle_formula,g0,tuple(psilist),tuple(kappalist),moduli_type)
  return answer

def socle_formula(g,psilist,kappalist,moduli_type=MODULI_ST):
  if moduli_type == MODULI_CT or g == 0:
    return CTconst(g)*CTsum(list(psilist),list(kappalist))
  if moduli_type <= MODULI_SM or moduli_type == MODULI_RT:
    return RTsum(g,list(psilist),list(kappalist))
  if moduli_type == MODULI_ST:
    return STsum(list(psilist),list(kappalist))

def multi(sigma):
  term = factorial(sum(sigma))
  for i in sigma:
    term /= factorial(i)
  return term

def multi2(g, sigma):
  sigma.sort()
  if sigma[0] == 0:
    total = 0
    for i in range(len(sigma)-1):
      sigmacopy = sigma[1:]
      if sigmacopy[i] > 0:
        sigmacopy[i] -= 1
        total += multi2(g,sigmacopy)
    return total
  term = factorial(2*g-3+len(sigma))
  term *= (2*g-1).multifactorial(2)
  term /= factorial(2*g-1)
  for i in sigma:
    term /= (2*i-1).multifactorial(2)
  return term

def STsum(psilist,kappalist):
  kappalist.sort()
  total = 0
  for i in setparts_with_auts(kappalist):
    total += (-1)^(len(i[0])) * i[1] * STrecur([1+sum(j) for j in i[0]] + psilist)
  return total*(-1)^(len(kappalist))

def RTsum(g,psilist,kappalist):
  kappalist.sort()
  total = 0
  for i in setparts_with_auts(kappalist):
    total += (-1)^(len(i[0])) * i[1] * multi2(g, [1+sum(j) for j in i[0]] + psilist)
  return total*(-1)^(len(kappalist))

def CTsum(psilist,kappalist):
  kappalist.sort()
  total = 0
  for i in setparts_with_auts(kappalist):
    total += (-1)^(len(i[0])) * i[1] * multi([1+sum(j) for j in i[0]] + psilist)
  return total*(-1)^(len(kappalist))

def CTconst(g):
  return abs((2^(2*g-1)-1)*bernoulli(2*g))/(2^(2*g-1)*factorial(2*g))

def STrecur(psi):
  psi.sort()
  return capply(STrecur_calc,tuple(psi))

def STrecur_calc(psi):
  n = len(psi)
  if n == 0:
    return 1
  s = sum(psi)
  if (s - n) % 3 != 0:
    return 0
  if psi[0] == 0:
    if s == 0 and n == 3:
      return 1
    total = 0
    for i in range(n-1):
      psicopy = list(psi[1:])
      if psicopy[i] > 0:
        psicopy[i] -= 1
        total += STrecur(psicopy)
    return total

  g = (s - n)/3 + 1
  d = psi[-1]
  total = 0

  psicopy = [0,0,0,0] + list(psi)
  psicopy[-1] += 1
  total += (2*d+3)/12 * STrecur(psicopy)

  psicopy = [0,0,0] + list(psi)
  total -= (2*g+n-1)/6 * STrecur(psicopy)

  for I in Subsets(range(n-1)):
    psi3 = [0,0] + [psi[i] for i in I]
    x = STrecur(psi3)
    if x == 0:
      continue
    psi1 = [0,0] + [psi[i] for i in range(n-1) if i not in I] + [d+1]
    psi2 = list(psi1[1:])
    psi2[-1] = d
    total += ((2*d+3)*STrecur(psi1) - (2*g+n-1)*STrecur(psi2))*x
  total /= (2*g+n-1)*(2*g+n-2)
  return total

def good_generator_list(g,r,markings=(),moduli_type=MODULI_ST):
  gens = capply(all_strata,g,r,markings,moduli_type)
  good_gens = []
  ngens = len(gens)
  for num in range(ngens):
    G = gens[num]
    good = True
    for i in range(1,G.M.nrows()):
      g = G.M[i,0][0]
      codim = 0
      for d in range(1,r+1):
        if G.M[i,0][d] != 0:
          if 3*d > g:
            good = False
            break
          codim += d*G.M[i,0][d]
      if not good:
        break
      for j in range(1,G.M.ncols()):
        codim += G.M[i,j][1]
        codim += G.M[i,j][2]
      if codim > 0 and codim >= g:
        good = False
        break
    if good:
      good_gens.append(num)
  return good_gens

def automorphism_cosets(num,g,r,markings=(),moduli_type=MODULI_ST):
  G = single_stratum(num,g,r,markings,moduli_type)
  pureG = Graph(G.M)
  pureG.purify()
  pureG.compute_invariant()
  pure_auts = graph_list_isomorphisms(pureG, pureG)
  num_pure = len(pure_auts)
  impure_auts = graph_list_isomorphisms(G, G)
  num_impure = len(impure_auts)
  chosen_auts = []
  used_auts = []
  v = G.num_vertices()
  e = G.num_edges()
  for i in range(num_pure):
    if i not in used_auts:
      chosen_auts.append(pure_auts[i])
      for g in impure_auts:
        sigma = [[pure_auts[i][0][g[0][k]-1] for k in range(v)],[pure_auts[i][1][g[1][k]-1] for k in range(e)]]
        for ii in range(num_pure):
          if pure_auts[ii] == sigma:
            used_auts.append(ii)
            break
  return [num_impure,chosen_auts]

def gorenstein(g,r,marked_points=(),moduli_type=MODULI_ST):
  """
  This function returns the rank of the codimension r grading of the
  Gorenstein quotient of the tautological ring of the moduli space of genus g
  curves with marked points labeled by the multiset marked_points.

  g and r should be nonnegative integers and marked_points should be a
  tuple of positive integers.

  The parameter moduli_type determines which moduli space to use:
  - MODULI_ST: all stable curves (this is the default)
  - MODULI_CT: curves of compact type
  - MODULI_RT: curves with rational tails

  EXAMPLES:

  - rank Gor^3(bar{M}_{3}) = 10
      sage: gorenstein(3,3)
      10

  - rank Gor^2(bar{M}_{2,2}) = 14
      sage: gorenstein(2,2,(1,2))           
      14

  - rank Gor^2(bar{M}_{2,2})^{S_2} = 11
      sage: gorenstein(2,2,(1,1))  
      11

  - rank Gor^2(M^c_{4}) = 6
      sage: gorenstein(4,2,(),MODULI_CT)
      6

  - rank Gor^4(M^rt_{8,2}) = 22
      sage: gorenstein(8,4,(1,2),MODULI_RT)
      22
  """
  #print "Computing lookup tables"
  gorenstein_precompute(g,r,marked_points,moduli_type)
  r3 = dim_form(g,len(marked_points),moduli_type)
  r2 = r3-r
  S1 = capply(good_generator_list,g,r,marked_points,moduli_type)
  S2 = capply(good_generator_list,g,r2,marked_points,moduli_type)
  M = pairing_submatrix(S1,S2,g,r,marked_points,moduli_type)
  #print "Computing rank"
  return matrix(M).rank()

def goren_rels(g,r,markings=(),moduli_type=MODULI_ST):
  gorenstein_precompute(g,r,markings,moduli_type)
  r3 = dim_form(g,len(markings),moduli_type)
  r2 = r3-r
  S1 = range(num_strata(g,r,markings,moduli_type))
  S2 = capply(good_generator_list,g,r2,markings,moduli_type)
  M = pairing_submatrix(S1,S2,g,r,markings,moduli_type)
  return Matrix(M).kernel().basis()

##### some crude parallelization:

@parallel
def para_pairing_dict(S1,S2,g,r1,markings=(),moduli_type=MODULI_ST):
  D = {}
  r3 = dim_form(g,len(markings),moduli_type)
  r2 = r3-r1
  for i1 in S1:
    for i2 in S2:
      D[i1,i2] = compute_pairing(r1,i1,r2,i2,g,r3,markings,moduli_type)
  return D

def para_gorenstein(g,r1,markings=(),moduli_type=MODULI_ST):
  r3 = dim_form(g,len(markings),moduli_type)
  r2 = r3-r1
  # precompute stuff
  D = capply(all_strata,g,r1,markings,moduli_type)
  D = capply(all_strata,g,r2,markings,moduli_type)
  D = capply(all_strata,g,r3,markings,moduli_type)
  for r in range(r3+1):
    D = capply(all_pure_strata,g,r,markings,moduli_type)
    for num in range(num_pure_strata(g,r,markings,moduli_type)):
      D = capply(pure_strata_autom_count,num,g,r,markings,moduli_type)
  D = capply(contraction_table,g,r3,markings,moduli_type)
  ngens3 = num_strata(g,r3,markings,moduli_type)
  for i in range(ngens3):
    D = capply(socle_evaluation,i,g,markings,moduli_type)
  # probably should do more stuff here

  # choose good generators
  S1 = capply(good_generator_list,g,r1,markings,moduli_type)
  S2 = capply(good_generator_list,g,r2,markings,moduli_type)

  # more precomputing
  for i1 in S1:
    D = capply(automorphism_cosets,i1,g,r1,markings,moduli_type)
  for i2 in S2:
    D = capply(automorphism_cosets,i2,g,r2,markings,moduli_type)

  # chop generator lists into blocks
  a1 = len(S1)
  a2 = len(S2)
  k1 = min(a1,10)
  k2 = min(a2,10)
  S1list = [S1[floor(i*a1/k1):floor((i+1)*a1/k1)] for i in range(k1)]
  S2list = [S2[floor(i*a2/k2):floor((i+1)*a2/k2)] for i in range(k2)]
  input_list = []
  for T1 in S1list:
    for T2 in S2list:
      input_list.append((T1,T2,g,r1,markings,moduli_type))
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