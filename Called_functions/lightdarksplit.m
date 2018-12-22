%takes in a vector and outputs all indices for daytime values and nighttime
%values
%based on 9AM-11PM light ON
%excludes first 2 hours of each phase

function [dayIndices, nightIndices] = lightdarksplit(timeVector);

nightIndices = intersect(find(timeVector>1),find(timeVector<9));
dayIndices = intersect(find(timeVector>11),find(timeVector<23));

end