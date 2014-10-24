def all_tests():
  clear_cache()
  strata_tests()
  betti_tests()
  gorenstein_tests()
  primrel_tests()
  recursive_betti_tests()

def test(f,*args):
  start_time = time.time()
  start_memory = floor(get_memory_usage())
  ans = apply(f,args)
  end_time = time.time()
  end_memory = floor(get_memory_usage())
  print "%s (T: %s, M: %s)" % (ans,floor(end_time-start_time),end_memory-start_memory)

def strata_tests():
  test(num_strata,5,6)
  test(num_strata,3,3,(1,2,3,4))
  test(num_strata,2,3,(1,1,1,1,1,1))
  test(num_strata,5,5,(1,2),MODULI_CT)
  test(num_strata,10,7,(1,1,1,1,1),MODULI_RT)

def betti_tests():
  test(betti,5,4)
  test(betti,2,3,(1,1,1))
  test(betti,1,2,(1,2,3,4))
  test(betti,5,5,(),MODULI_CT)
  test(betti,16,8,(1,),MODULI_SM)

def gorenstein_tests():
  test(gorenstein,3,2,(1,))
  test(gorenstein,3,3,(1,))
  test(gorenstein,5,3,(),MODULI_CT)
  test(gorenstein,2,2,(1,1,1,1),MODULI_CT)
  test(gorenstein,16,7,(1,),MODULI_RT)

def primrel_tests():
  test(num_new_rels,6,4)
  test(num_new_rels,2,3,4)
  test(num_new_rels,3,3,3)
  test(num_new_rels,4,3,3)
  test(num_new_rels,13,7,0,MODULI_SM)

def recursive_betti_tests():
  test(recursive_betti,5,4)
  test(recursive_betti,2,3,(1,1,1))
  test(recursive_betti,1,2,(1,2,3,4))
  test(recursive_betti,5,5,(),MODULI_CT)
  test(recursive_betti,16,8,(1,),MODULI_SM)