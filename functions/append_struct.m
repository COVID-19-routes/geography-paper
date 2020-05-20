function S1=append_struct(S1,S2)
%append fields of structure S2 to structure S1;

f = fieldnames(S2);
for i = 1:length(f)
    S1.(f{i})=S2.(f{i});
end