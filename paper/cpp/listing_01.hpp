#define return_macro(expr) \
  -> decltype(safeToReturn(expr)) \
  { return safeToReturn(expr); } 
