% Overloaded Division Operator for Parameters Structure

function c = mrdivide(p,q)
ui = 'SimBiology';
% Check for non-structs
if (isa(p,'numeric'))
    a.Value = p;
    a.Units = 'dimensionless';
else
    a = p;
end
if (isa(q,'numeric'))
    b.Value = q;
    b.Units = 'dimensionless';
else
    b = q;
end
c.Value = a.Value/b.Value;
c.Units = symunit2str(str2symunit(a.Units,ui)/str2symunit(b.Units,ui),ui);
c = simplify(c);
end