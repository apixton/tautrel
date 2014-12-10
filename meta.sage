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

def dsave(str,*args):
  if ENABLE_DSAVE:
    save(0,str % args)