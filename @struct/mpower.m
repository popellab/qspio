% Overloaded Power Operator for Parameters Structure

function c = mpower(a,n)
ui = 'SimBiology';
c.Value = a.Value^n;
c.Units = symunit2str(str2symunit(a.Units,ui)^n,ui);
end