def application_triqs_import(name,*args,**kwargs):
  if name.startswith('@package_name@'):
    name = name[len('@package_name@')+1:]
  return builtin_import(name,*args,**kwargs)

import builtins
builtins.__import__, builtin_import = application_triqs_import, builtins.__import__

