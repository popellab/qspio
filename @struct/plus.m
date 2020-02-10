% Overloaded Multiplication for Parameter Structures

function c = plus(a,b)

% Get Symbolic Units
ua = str2symunit(a.Units,'SimBiology');
ub = str2symunit(b.Units,'SimBiology');

% Get Conversion Factor
fac = double(unitConversionFactor(ua,ub));

% Perform Addition
c.Value = a.Value + b.Value/fac;

% Tranfer Units
c.Units = a.Units;