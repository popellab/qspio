% Overloaded 'simplify' for Parameter Structures

function b = simflify(a)
  
  ui = 'SimBiology';
  
  u = str2symunit(a.Units,ui);
  expr = simplify(a.Value*u);
  [value,units] =  separateUnits(expr);
  b.Value = double(value);
  b.Units = symunit2str(units,ui);
  
end