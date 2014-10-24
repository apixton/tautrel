def all_tests():
  clear_cache()
  strata_tests()
  betti_tests()
  gorenstein_tests()
  primrel_tests()
  recursive_betti_tests()

def strata_tests():
  print num_strata(5,6)
  print num_strata(3,3,(1,2,3,4))
  print num_strata(2,3,(1,1,1,1,1,1))
  print num_strata(5,5,(1,2),MODULI_CT)
  print num_strata(10,7,(1,1,1,1,1),MODULI_RT)

def betti_tests():
  print betti(5,4)
  print betti(2,3,(1,1,1))
  print betti(1,2,(1,2,3,4))
  print betti(5,5,(),MODULI_CT)
  print betti(16,8,(1,),MODULI_SM)

def gorenstein_tests():
  print gorenstein(3,2,(1,))
  print gorenstein(3,3,(1,))
  print gorenstein(5,3,(),MODULI_CT)
  print gorenstein(2,2,(1,1,1,1),MODULI_CT)
  print gorenstein(16,7,(1,),MODULI_RT)

def primrel_tests():
  print num_new_rels(6,4)
  print num_new_rels(2,3,4)
  print num_new_rels(3,3,3)
  print num_new_rels(4,3,3)
  print num_new_rels(13,7,0,MODULI_SM)

def recursive_betti_tests():
  print recursive_betti(5,4)
  print recursive_betti(2,3,(1,1,1))
  print recursive_betti(1,2,(1,2,3,4))
  print recursive_betti(5,5,(),MODULI_CT)
  print recursive_betti(16,8,(1,),MODULI_SM)