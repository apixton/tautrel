import time
import os.path

def clear_cache():
  global cache_dict
  cache_dict = {}

def capply(f,*args):
  global cache_dict
  if not cache_dict.has_key(f.func_name):
    cache_dict[f.func_name] = {}
  elif cache_dict[f.func_name].has_key(args):
    return cache_dict[f.func_name][args]
  fc_list = ['choose_basic_rels','all_strata','partial_symmetrize_map',
             #'single_insertion_pullback','single_insertion_pullback2',
             #'single_kappa_multiple','single_psi_multiple',
             'strata_invariant_lookup']
  if f.func_name in fc_list:
    return fcapply(f,*args)
  cache_dict[f.func_name][args] = apply(f,args)
  return cache_dict[f.func_name][args]

def fcapply(f,*args):
  global cache_dict
  if not cache_dict.has_key(f.func_name):
    cache_dict[f.func_name] = {}
  elif cache_dict[f.func_name].has_key(args):
    return cache_dict[f.func_name][args]
  filename = 'obj/' + f.func_name + ':' + ':'.join(fixup_args(f,args))
  if os.path.isfile(filename + '.sobj'):
    ans = load(filename)
    cache_dict[f.func_name][args] = ans
    return ans
  ans = apply(f,args)
  cache_dict[f.func_name][args] = ans
  if not os.path.isfile(filename + '.sobj'):
    save(ans,filename)
  return ans

def dlog(files,str,*args):
  if isinstance(files,list):
    filelist = files
  else:
    filelist = [files]
  for filename in filelist:
    o = open('logs/' + filename + '.txt', 'a')
    o.write(time.ctime() + ': ' + (str % args) + '\n')
    o.close()

def mod_type_string(moduli_type):
  if moduli_type == MODULI_ST:
    return "MODULI_ST"
  elif moduli_type == MODULI_CT:
    return "MODULI_CT"
  elif moduli_type == MODULI_RT:
    return "MODULI_RT"
  elif moduli_type == MODULI_SM:
    return "MODULI_SM"
  return "MODULI_XX"

def fixup_args(f,arg_list):
  new_list = []
  count = 0
  for arg in arg_list:
    if f.func_code.co_varnames[count] == 'moduli_type':
      new_list.append(mod_type_string(arg))
    else:
      new_list.append(str(arg))
    count += 1
  return new_list

def log_func(f,*args):
  comp = f.func_name + '(' + ','.join(fixup_args(f,args)) + ')'
  dlog('history','began computing %s',comp)
  start_time = time.time()
  start_memory = floor(get_memory_usage())
  result = apply(f,args)
  end_time = time.time()
  end_memory = floor(get_memory_usage())
  dlog(['results','history'],'%s = %s',comp,result)
  dlog('history','finished computing %s: %s sec, %s MB',comp,floor(end_time-start_time),end_memory-start_memory)
  return result

def random_permutation(L):
  l = len(L)
  sigma = Permutations(l).random_element()
  LL = [L[i-1] for i in sigma]
  return LL