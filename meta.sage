import time

def clear_cache():
  global cache_dict
  cache_dict = {}

def capply(f,*args):
  global cache_dict
  if not cache_dict.has_key(f.func_name):
    cache_dict[f.func_name] = {}
  if not cache_dict[f.func_name].has_key(args):
    cache_dict[f.func_name][args] = apply(f,args)
  return cache_dict[f.func_name][args]

def dprint(str,*args):
  if ENABLE_DPRINT:
    print str % args

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