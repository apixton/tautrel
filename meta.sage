def clear_cache():
  global cache_dict
  cache_dict = {}

def capply(f,*args):
  global cache_dict
  key = (f.func_name,args)
  if not cache_dict.has_key(key):
    cache_dict[key] = apply(f,args)
  return cache_dict[key]

def dprint(str,*args):
  if ENABLE_DPRINT:
    print str % args

def dsave(str,*args):
  if ENABLE_DSAVE:
    save(0,str % args)